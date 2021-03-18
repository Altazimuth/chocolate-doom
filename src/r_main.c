// Emacs style mode select   -*- C++ -*-
//-----------------------------------------------------------------------------
//
// Copyright(C) 1993-1996 Id Software, Inc.
// Copyright(C) 2005 Simon Howard
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
// 02111-1307, USA.
//
// DESCRIPTION:
//  Rendering main loop and setup functions,
//   utility functions (BSP, geometry, trigonometry).
//  See tables.c, too.
//
//-----------------------------------------------------------------------------





#include <stdlib.h>
#include <math.h>


#include "doomdef.h"
#include "d_net.h"

#include "m_bbox.h"
#include "m_menu.h"

#include "r_local.h"
#include "r_sky.h"

#include "v_video.h"
#include "doomstat.h"



// Fineangles in the SCREENWIDTH wide window.
#define FIELDOFVIEW		2048



int viewangleoffset;

// increment every time a check is made
int validcount = 1;


lighttable_t* fixedcolormap;
extern lighttable_t** walllights;

int centerx;
int centery;

fixed_t centerxfrac;
fixed_t centeryfrac;
fixed_t projection;

// just for profiling purposes
int framecount;

int sscount;
int linecount;
int loopcount;

fixed_t viewx;
fixed_t viewy;
fixed_t viewz;

angle_t viewangle;

fixed_t viewcos;
fixed_t viewsin;

player_t* viewplayer;

// 0 = high, 1 = low
int detailshift;

//
// precalculated math tables
//
angle_t clipangle;

// The viewangletox[viewangle + FINEANGLES/4] lookup
// maps the visible view angles to screen X coordinates,
// flattening the arc to a flat projection plane.
// There will be many angles mapped to the same X.
int viewangletox[FINEANGLES / 2];

// The xtoviewangleangle[] table maps a screen pixel
// to the lowest viewangle that maps back to x ranges
// from clipangle to -clipangle.
angle_t xtoviewangle[SCREENWIDTH + 1];


// UNUSED.
// The finetangentgent[angle+FINEANGLES/4] table
// holds the fixed_t tangent values for view angles,
// ranging from INT_MIN to 0 to INT_MAX.
// fixed_t      finetangent[FINEANGLES/2];

// fixed_t      finesine[5*FINEANGLES/4];
const fixed_t* finecosine = &finesine[FINEANGLES / 4];


lighttable_t* scalelight[LIGHTLEVELS][MAXLIGHTSCALE];
lighttable_t* scalelightfixed[MAXLIGHTSCALE];
lighttable_t* zlight[LIGHTLEVELS][MAXLIGHTZ];

// bumped light from gun blasts
int extralight;



void (*colfunc) (void);
void (*basecolfunc) (void);
void (*fuzzcolfunc) (void);
void (*transcolfunc) (void);
void (*spanfunc) (visplane_t* pl, visplane_t* firstplane);



//
// R_AddPointToBox
// Expand a given bbox
// so that it encloses a given point.
//
void R_AddPointToBox(int x, int y, fixed_t* box)
{
	if (x < box[BOXLEFT])
		box[BOXLEFT] = x;
	if (x > box[BOXRIGHT])
		box[BOXRIGHT] = x;
	if (y < box[BOXBOTTOM])
		box[BOXBOTTOM] = y;
	if (y > box[BOXTOP])
		box[BOXTOP] = y;
}


//
// R_PointOnSide
// Traverse BSP (sub) tree,
//  check point against partition plane.
// Returns side 0 (front) or 1 (back).
//
int R_PointOnSide(fixed_t x, fixed_t y, node_t* node)
{
	fixed_t dx;
	fixed_t dy;
	fixed_t left;
	fixed_t right;
	
	if (!node->dx)
	{
		if (x <= node->x)
			return node->dy > 0;
			
		return node->dy < 0;
	}
	if (!node->dy)
	{
		if (y <= node->y)
			return node->dx < 0;
			
		return node->dx > 0;
	}
	
	dx = (x - node->x);
	dy = (y - node->y);
	
	// Try to quickly decide by looking at sign bits.
	if ((node->dy ^ node->dx ^ dx ^ dy) & 0x80000000)
	{
		if ((node->dy ^ dx) & 0x80000000)
		{
			// (left is negative)
			return 1;
		}
		return 0;
	}
	
	left = FixedMul(node->dy >> FRACBITS, dx);
	right = FixedMul(dy, node->dx >> FRACBITS);
	
	if (right < left)
	{
		// front side
		return 0;
	}
	// back side
	return 1;
}


int R_PointOnSegSide(fixed_t x, fixed_t y, seg_t* line)
{
	fixed_t lx;
	fixed_t ly;
	fixed_t ldx;
	fixed_t ldy;
	fixed_t dx;
	fixed_t dy;
	fixed_t left;
	fixed_t right;
	
	lx = line->v1->x;
	ly = line->v1->y;
	
	ldx = line->v2->x - lx;
	ldy = line->v2->y - ly;
	
	if (!ldx)
	{
		if (x <= lx)
			return ldy > 0;
			
		return ldy < 0;
	}
	if (!ldy)
	{
		if (y <= ly)
			return ldx < 0;
			
		return ldx > 0;
	}
	
	dx = (x - lx);
	dy = (y - ly);
	
	// Try to quickly decide by looking at sign bits.
	if ((ldy ^ ldx ^ dx ^ dy) & 0x80000000)
	{
		if ((ldy ^ dx) & 0x80000000)
		{
			// (left is negative)
			return 1;
		}
		return 0;
	}
	
	left = FixedMul(ldy >> FRACBITS, dx);
	right = FixedMul(dy, ldx >> FRACBITS);
	
	if (right < left)
	{
		// front side
		return 0;
	}
	// back side
	return 1;
}


//
// R_PointToAngle
// To get a global angle from cartesian coordinates,
//  the coordinates are flipped until they are in
//  the first octant of the coordinate system, then
//  the y (<=x) is scaled and divided by x to get a
//  tangent (slope) value which is looked up in the
//  tantoangle[] table.

//




angle_t R_PointToAngle(fixed_t x, fixed_t y)
{
	x -= viewx;
	y -= viewy;
	
	if ((!x) && (!y))
		return 0;
		
	if (x >= 0)
	{
		// x >=0
		if (y >= 0)
		{
			// y>= 0
			
			if (x > y)
			{
				// octant 0
				return tantoangle[SlopeDiv(y, x)];
			}
			else
			{
				// octant 1
				return ANG90 - 1 - tantoangle[SlopeDiv(x, y)];
			}
		}
		else
		{
			// y<0
			y = -y;
			
			if (x > y)
			{
				// octant 8
				return -tantoangle[SlopeDiv(y, x)];
			}
			else
			{
				// octant 7
				return ANG270 + tantoangle[SlopeDiv(x, y)];
			}
		}
	}
	else
	{
		// x<0
		x = -x;
		
		if (y >= 0)
		{
			// y>= 0
			if (x > y)
			{
				// octant 3
				return ANG180 - 1 - tantoangle[SlopeDiv(y, x)];
			}
			else
			{
				// octant 2
				return ANG90 + tantoangle[SlopeDiv(x, y)];
			}
		}
		else
		{
			// y<0
			y = -y;
			
			if (x > y)
			{
				// octant 4
				return ANG180 + tantoangle[SlopeDiv(y, x)];
			}
			else
			{
				// octant 5
				return ANG270 - 1 - tantoangle[SlopeDiv(x, y)];
			}
		}
	}
	return 0;
}


angle_t R_PointToAngle2(fixed_t x1, fixed_t y1, fixed_t x2, fixed_t y2)
{
	viewx = x1;
	viewy = y1;
	
	return R_PointToAngle(x2, y2);
}


fixed_t R_PointToDist(fixed_t x, fixed_t y)
{
	int angle;
	fixed_t dx;
	fixed_t dy;
	fixed_t temp;
	fixed_t dist;
	fixed_t frac;
	
	dx = abs(x - viewx);
	dy = abs(y - viewy);
	
	if (dy > dx)
	{
		temp = dx;
		dx = dy;
		dy = temp;
	}
	// Fix crashes in udm1.wad
	
	if (dx != 0)
	{
		frac = FixedDiv(dy, dx);
	}
	else
	{
		frac = 0;
	}
	
	angle = (tantoangle[frac >> DBITS] + ANG90) >> ANGLETOFINESHIFT;
	
	// use as cosine
	dist = FixedDiv(dx, finesine[angle]);
	
	return dist;
}




//
// R_InitPointToAngle
//
void R_InitPointToAngle(void)
{
	// UNUSED - now getting from tables.c
#if 0
	int i;
	long t;
	float f;
	
//
// slope (tangent) to angle lookup
//
	for (i = 0; i <= SLOPERANGE; i++)
	{
		f = atan((float)i / SLOPERANGE) / (3.141592657 * 2);
		t = 0xffffffff * f;
		tantoangle[i] = t;
	}
#endif
}


//
// R_ScaleFromGlobalAngle
// Returns the texture mapping scale
//  for the current line (horizontal span)
//  at the given angle.
// rw_distance must be calculated first.
//
fixed_t R_ScaleFromGlobalAngle(angle_t visangle)
{
	fixed_t scale;
	angle_t anglea;
	angle_t angleb;
	int sinea;
	int sineb;
	fixed_t num;
	int den;
	
	// UNUSED
#if 0
	{
		fixed_t dist;
		fixed_t z;
		fixed_t sinv;
		fixed_t cosv;
		
		sinv = finesine[(visangle - rw_normalangle) >> ANGLETOFINESHIFT];
		dist = FixedDiv(rw_distance, sinv);
		cosv = finecosine[(viewangle - visangle) >> ANGLETOFINESHIFT];
		z = abs(FixedMul(dist, cosv));
		scale = FixedDiv(projection, z);
		return scale;
	}
#endif
	
	anglea = ANG90 + (visangle - viewangle);
	angleb = ANG90 + (visangle - rw_normalangle);
	
	// both sines are allways positive
	sinea = finesine[anglea >> ANGLETOFINESHIFT];
	sineb = finesine[angleb >> ANGLETOFINESHIFT];
	num = FixedMul(projection, sineb) << detailshift;
	den = FixedMul(rw_distance, sinea);
	
	if (den > num >> 16)
	{
		scale = FixedDiv(num, den);
		
		if (scale > 64 * FRACUNIT)
			scale = 64 * FRACUNIT;
		else if (scale < 256)
			scale = 256;
	}
	else
		scale = 64 * FRACUNIT;
		
	return scale;
}



//
// R_InitTables
//
void R_InitTables(void)
{
	// UNUSED: now getting from tables.c
#if 0
	int i;
	float a;
	float fv;
	int t;
	
	// viewangle tangent table
	for (i = 0; i < FINEANGLES / 2; i++)
	{
		a = (i - FINEANGLES / 4 + 0.5) * PI * 2 / FINEANGLES;
		fv = FRACUNIT * tan(a);
		t = fv;
		finetangent[i] = t;
	}
	
	// finesine table
	for (i = 0; i < 5 * FINEANGLES / 4; i++)
	{
		// OPTIMIZE: mirror...
		a = (i + 0.5) * PI * 2 / FINEANGLES;
		t = FRACUNIT * sin(a);
		finesine[i] = t;
	}
#endif
	
}



//
// R_InitTextureMapping
//
void R_InitTextureMapping(void)
{
	int i;
	int x;
	int t;
	fixed_t focallength;
	
	// Use tangent table to generate viewangletox:
	//  viewangletox will give the next greatest x
	//  after the view angle.
	//
	// Calc focallength
	//  so FIELDOFVIEW angles covers SCREENWIDTH.
	focallength = FixedDiv(centerxfrac, finetangent[FINEANGLES / 4 + FIELDOFVIEW / 2]);
	
	for (i = 0; i < FINEANGLES / 2; i++)
	{
		if (finetangent[i] > FRACUNIT * 2)
			t = -1;
		else if (finetangent[i] < -FRACUNIT * 2)
			t = viewwidth + 1;
		else
		{
			t = FixedMul(finetangent[i], focallength);
			t = (centerxfrac - t + FRACUNIT - 1) >> FRACBITS;
			
			if (t < -1)
				t = -1;
			else if (t > viewwidth + 1)
				t = viewwidth + 1;
		}
		viewangletox[i] = t;
	}
	
	// Scan viewangletox[] to generate xtoviewangle[]:
	//  xtoviewangle will give the smallest view angle
	//  that maps to x.
	for (x = 0; x <= viewwidth; x++)
	{
		i = 0;
		while (viewangletox[i] > x)
			i++;
		xtoviewangle[x] = (i << ANGLETOFINESHIFT) - ANG90;
	}
	
	// Take out the fencepost cases from viewangletox.
	for (i = 0; i < FINEANGLES / 2; i++)
	{
		t = FixedMul(finetangent[i], focallength);
		t = centerx - t;
		
		if (viewangletox[i] == -1)
			viewangletox[i] = 0;
		else if (viewangletox[i] == viewwidth + 1)
			viewangletox[i] = viewwidth;
	}
	
	clipangle = xtoviewangle[0];
}



//
// R_InitLightTables
// Only inits the zlight table,
//  because the scalelight table changes with view size.
//
#define DISTMAP		2

void R_InitLightTables(void)
{
	int i;
	int j;
	int level;
	int startmap;
	int scale;
	
	// Calculate the light levels to use
	//  for each level / distance combination.
	for (i = 0; i < LIGHTLEVELS; i++)
	{
		startmap = ((LIGHTLEVELS - 1 - i) * 2) * NUMCOLORMAPS / LIGHTLEVELS;
		for (j = 0; j < MAXLIGHTZ; j++)
		{
			scale = FixedDiv((SCREENWIDTH / 2 * FRACUNIT), (j + 1) << LIGHTZSHIFT);
			scale >>= LIGHTSCALESHIFT;
			level = startmap - scale / DISTMAP;
			
			if (level < 0)
				level = 0;
				
			if (level >= NUMCOLORMAPS)
				level = NUMCOLORMAPS - 1;
				
			zlight[i][j] = colormaps + level * 256;
		}
	}
}



//
// R_SetViewSize
// Do not really change anything here,
//  because it might be in the middle of a refresh.
// The change will take effect next refresh.
//
boolean setsizeneeded;
int setblocks;
int setdetail;


void R_SetViewSize(int blocks, int detail)
{
	setsizeneeded = true;
	setblocks = blocks;
	setdetail = detail;
}


//
// R_ExecuteSetViewSize
//
void R_ExecuteSetViewSize(void)
{
	fixed_t cosadj;
	fixed_t dy;
	int i;
	int j;
	int level;
	int startmap;
	
	setsizeneeded = false;
	
	if (setblocks == 11)
	{
		scaledviewwidth = SCREENWIDTH;
		viewheight = SCREENHEIGHT;
	}
	else
	{
		scaledviewwidth = setblocks * 32;
		viewheight = (setblocks * 168 / 10) & ~7;
	}
	
	detailshift = setdetail;
	viewwidth = scaledviewwidth >> detailshift;
	
	centery = viewheight / 2;
	centerx = viewwidth / 2;
	centerxfrac = centerx << FRACBITS;
	centeryfrac = centery << FRACBITS;
	projection = centerxfrac;
	
	if (!detailshift)
	{
		colfunc = basecolfunc = R_DrawColumn;
		fuzzcolfunc = R_DrawFuzzColumn;
		transcolfunc = R_DrawTranslatedColumn;
		spanfunc = R_DrawSpan;
	}
	else
	{
		colfunc = basecolfunc = R_DrawColumnLow;
		fuzzcolfunc = R_DrawFuzzColumnLow;
		transcolfunc = R_DrawTranslatedColumnLow;
		spanfunc = R_DrawSpanLow;
	}
	
	R_InitBuffer(scaledviewwidth, viewheight);
	
	R_InitTextureMapping();
	
	// psprite scales
	pspritescale = FRACUNIT * viewwidth / SCREENWIDTH;
	pspriteiscale = FRACUNIT * SCREENWIDTH / viewwidth;
	
	// thing clipping
	for (i = 0; i < viewwidth; i++)
		screenheightarray[i] = viewheight;
		
	// planes
	for (i = 0; i < viewheight; i++)
	{
		dy = ((i - viewheight / 2) << FRACBITS) + FRACUNIT / 2;
		dy = abs(dy);
		yslope[i] = FixedDiv((viewwidth << detailshift) / 2 * FRACUNIT, dy);
	}
	
	for (i = 0; i < viewwidth; i++)
	{
		cosadj = abs(finecosine[xtoviewangle[i] >> ANGLETOFINESHIFT]);
		distscale[i] = FixedDiv(FRACUNIT, cosadj);
	}
	
	// Calculate the light levels to use
	//  for each level / scale combination.
	for (i = 0; i < LIGHTLEVELS; i++)
	{
		startmap = ((LIGHTLEVELS - 1 - i) * 2) * NUMCOLORMAPS / LIGHTLEVELS;
		for (j = 0; j < MAXLIGHTSCALE; j++)
		{
			level = startmap - j * SCREENWIDTH / (viewwidth << detailshift) / DISTMAP;
			
			if (level < 0)
				level = 0;
				
			if (level >= NUMCOLORMAPS)
				level = NUMCOLORMAPS - 1;
				
			scalelight[i][j] = colormaps + level * 256;
		}
	}
}



//
// R_Init
//



void R_Init(void)
{
	R_InitData();
	printf(".");
	R_InitPointToAngle();
	printf(".");
	R_InitTables();
	// viewwidth / viewheight / detailLevel are set by the defaults
	printf(".");
	
	R_SetViewSize(screenblocks, detailLevel);
	R_InitPlanes();
	printf(".");
	R_InitLightTables();
	printf(".");
	R_InitSkyMap();
	R_InitTranslationTables();
	printf(".");
	
	framecount = 0;
}


//
// R_PointInSubsector
//
subsector_t* R_PointInSubsector(fixed_t x, fixed_t y)
{
	node_t* node;
	int side;
	int nodenum;
	
	// single subsector is a special case
	if (!numnodes)
		return subsectors;
		
	nodenum = numnodes - 1;
	
	while (!(nodenum & NF_SUBSECTOR))
	{
		node = &nodes[nodenum];
		side = R_PointOnSide(x, y, node);
		nodenum = node->children[side];
	}
	
	return &subsectors[nodenum & ~NF_SUBSECTOR];
}


//
// R_IsPointInSubsector, same of above but return 0 if not in subsector
//
subsector_t* R_IsPointInSubsector(fixed_t x, fixed_t y)
{
	node_t* node;
	int side;
	int nodenum, i;
	subsector_t* ret;
	
	// single subsector is a special case
	if (!numnodes)
		return subsectors;
		
	nodenum = numnodes - 1;
	
	while (!(nodenum & NF_SUBSECTOR))
	{
		node = &nodes[nodenum];
		side = R_PointOnSide(x, y, node);
		nodenum = node->children[side];
	}
	
	ret = &subsectors[nodenum & ~NF_SUBSECTOR];
	for (i = 0; i < ret->numlines; i++)
	{
		if (R_PointOnSegSide(x, y, &segs[ret->firstline + i]))
			return 0;
	}
	
	return ret;
}


//
// R_SetupFrame
//
void R_SetupFrame(player_t* player)
{
	int i;
	
	viewplayer = player;
	viewx = player->mo->x;
	viewy = player->mo->y;
	viewangle = player->mo->angle + viewangleoffset;
	extralight = player->extralight;
	
	viewz = player->viewz;
	
	viewsin = finesine[viewangle >> ANGLETOFINESHIFT];
	viewcos = finecosine[viewangle >> ANGLETOFINESHIFT];
	
	sscount = 0;
	
	if (player->fixedcolormap)
	{
		fixedcolormap = colormaps + player->fixedcolormap * 256 * sizeof(lighttable_t);
		
		walllights = scalelightfixed;
		
		for (i = 0; i < MAXLIGHTSCALE; i++)
			scalelightfixed[i] = fixedcolormap;
	}
	else
		fixedcolormap = 0;
		
	framecount++;
	validcount++;
}

extern int g_NoFlashyPurpleHOM;

static void DrawHOMBox(void)
{
	int color;
	
	// Draw flashing box over the window where rendering takes place.
	// This will highlight any HOM areas.
	
	if (g_NoFlashyPurpleHOM)
		color = 0;
	else
		color = 250 + (framecount % 5);
	V_DrawRect(viewwindowx, viewwindowy, viewwidth, viewheight, color);
}

uint8_t PlaneBuffer[SCREENWIDTH* SCREENHEIGHT + 1];

//
// R_RenderView
//
void R_RenderPlayerView(player_t* player)
{
	int i;
	
	R_SetupFrame(player);
	DrawHOMBox();
	
	// Clear buffers.
	R_ClearClipSegs();
	R_ClearDrawSegs();
	R_ClearPlanes();
	R_ClearSprites();
	
	// Clear plane buffer
	memset(PlaneBuffer, 0xFF, sizeof(PlaneBuffer));
	
	// GhostlyDeath <September 23, 2011> -- Clear sneakies
	if (!g_crlsneakmode)
	{
		for (i = 0; i < numlines; i++)
			lines[i].Sneaky = 0;
		for (i = 0; i < numsides; i++)
			sides[i].Sneaky = 0;
		for (i = 0; i < numsectors; i++)
			sectors[i].Sneaky = sectors[i].Flashy = 0;
		for (i = 0; i < numsubsectors; i++)
			subsectors[i].Sneaky = subsectors[i].Flashy = 0;
	}
	// check for new console commands.
	NetUpdate();
	
	// The head node is the last node output.
	R_RenderBSPNode(numnodes - 1);
	
	// Check for new console commands.
	NetUpdate();
	
	R_DrawPlanes();
	
	// Check for new console commands.
	NetUpdate();
	
	R_DrawMasked();
	
	// GhostlyDeath <September 23, 2011> -- Draw planes
	R_RenderPlaneOverlay();
	
	// Check for new console commands.
	NetUpdate();
}

/* R_HeatStruct_t -- Heat structure */
typedef struct R_HeatStruct_s
{
	uint8_t NoWhere;							// Spot is nowhere!
	int16_t XPos, YPos;							// Position
	int32_t VisPlanes;							// Visplane count
	int32_t DrawSegs;							// Draw segs
	int32_t Openings;							// Openings
	int32_t Sprites;							// Sprites
} R_HeatStruct_t;

extern int ActualDrawSegs;

/* R_BruteForce() - Brute force sectors */
void R_BruteForce(
	const int a_FirstSec,			// -brutestart
	const int a_LastSec,			// -brutestop
	fixed_t FRACSTEP,				// -brutegran
	const char* const a_DATFile,	// -brutedatfile
	const char* const a_PPMFile,	// -bruteppmfile
	const char* const a_LogFile,	// -brutelogfile
	
	int32_t a_BruteStep,			// -brutestep
	int32_t a_BruteStepOffset,		// -brutestepoffset
	
	fixed_t a_bbX1,					// -brutebounds
	fixed_t a_bbY1,
	fixed_t a_bbX2,
	fixed_t a_bbY2
	)
{
	char Buf[128];
	int s, Ang, zzPos;
	angle_t RealAng;
	fixed_t x, y, z;
	fixed_t sx, sy;
	fixed_t ex, ey;
	fixed_t bsx, bsy;
	fixed_t bex, bey;
	fixed_t OldX, OldY, OldZ;
	sector_t* Sec;
	subsector_t* SubS, *RealSubS;
	size_t RenderRun;
	uint32_t u32, i, j;
	uint8_t u8r, u8g, u8b;
	int RenderTheView;
	int16_t sixx, sixy;
	double dx, dy, ddd, ddz;
	
	uint32_t PosChecks, SawPlanes, SawSegs, SawOpenings, SawSprites;
	
	fixed_t SuperBB[4];
	
	fixed_t tMx, tMy;
	int32_t ImW, ImH, SpotX, SpotY, Temp;
	R_HeatStruct_t* ImHeats;
	R_HeatStruct_t* ThisHeat;
	
	FILE* RawData, *NicePic;
	
	FILE* LogOut;
	
	/* Open Log File */
	if (!a_LogFile)
		LogOut = stdout;
	else
	{
		LogOut = fopen(a_LogFile, "wt");
		
		if (!LogOut)
			LogOut = stdout;
		else
		{
			// Make it flush more
			setvbuf(LogOut, (char*)NULL, _IONBF, 0);
		}
	}
	
	/* Default gran? */
	if (FRACSTEP <= 0)
		FRACSTEP = (32 << FRACBITS);
	
	/* Default Brute Step? */
	if (!a_BruteStep)
		a_BruteStep = 1;
	
	// Limit
	if (a_BruteStep < 1)
		a_BruteStep = 1;
	
	/* Unset player position */
	OldX = players[0].mo->x;
	OldY = players[0].mo->y;
	OldZ = players[0].mo->z;
	P_UnsetThingPosition(players[0].mo);
	
	/* First determine the bounds of everything */
	// Clear big box
	M_ClearBox(SuperBB);
	
	// Using pre-determined bounding box?
	if (a_bbX1 || a_bbY1 || a_bbX2 || a_bbY2)
	{
		// Just use those coordinates instead
		M_AddToBox(&SuperBB, a_bbX1, a_bbY1);
		M_AddToBox(&SuperBB, a_bbX2, a_bbY2);
	}
	
	// Use based on input and output sectors
	else
		for (s = a_FirstSec; s < a_LastSec; s++)
		{
			// Get sector
			Sec = &sectors[s];
		
			// Add bounds to big box
			M_AddToBox(&SuperBB, Sec->BBox[BOXLEFT], Sec->BBox[BOXTOP]);
			M_AddToBox(&SuperBB, Sec->BBox[BOXRIGHT], Sec->BBox[BOXBOTTOM]);
		}
	
	/* Round bounding box coordinates (and origin) */
	for (i = 0; i < 4; i++)
#if 1
	{
		ddz = FIXED_TO_FLOAT(FRACSTEP);
		
		ddd = FIXED_TO_FLOAT(SuperBB[i]);
		ddd = trunc(((double)((int)(ddd / ddz))) * ddz);
		SuperBB[i] = FLOAT_TO_FIXED(ddd);
	}
#else
		SuperBB[i] = FixedMul(((FixedDiv(SuperBB[i], FRACSTEP) >> FRACBITS) << FRACBITS), FRACSTEP);
#endif
	
	/* Allocation Loop */
	// Allocate
	do
	{
		// Determine image size
		ImW = abs(FixedDiv(SuperBB[BOXRIGHT] - SuperBB[BOXLEFT], FRACSTEP) >> FRACBITS) / a_BruteStep;
		ImH = abs(FixedDiv(SuperBB[BOXTOP] - SuperBB[BOXBOTTOM], FRACSTEP) >> FRACBITS);
	
		ImHeats = malloc(sizeof(*ImHeats) * (ImW * ImH));
		
		// Failed?
		if (!ImHeats)
		{
			fprintf(LogOut, "CRL: Not enough memory. Lowering granularity.\n");
			FRACSTEP = FixedMul(FRACSTEP, 2 << FRACBITS);
			
			// Really big?
			if (FRACSTEP > (8192 << FRACBITS))
			{
				fprintf(LogOut, "CRL: Granularity too low.\n");
				return;
			}
		}
		
		fflush(LogOut);
	} while (!ImHeats);

	memset(ImHeats, 0, sizeof(*ImHeats) * (ImW * ImH));
	
	// Determine map origin
	tMx = SuperBB[BOXRIGHT];
	if (SuperBB[BOXLEFT] < tMx)
		tMx = SuperBB[BOXLEFT];
		
	tMy = SuperBB[BOXBOTTOM];
	if (SuperBB[BOXTOP] < tMy)
		tMy = SuperBB[BOXTOP];
		
	/* Determine bounds to check */
	bsx = SuperBB[BOXLEFT];
	bex = SuperBB[BOXRIGHT];
	bsy = SuperBB[BOXBOTTOM];
	bey = SuperBB[BOXTOP];
	
	// Correct?
	if (bex < bsx)
	{
		z = bex;
		bex = bsx;
		bsx = z;
	}
	
	if (bey < bsy)
	{
		z = bey;
		bey = bsy;
		bsy = z;
	}
	
	/* Go through all sectors */
	RenderRun = PosChecks = SawPlanes = SawSegs = SawOpenings = SawSprites = 0;
	for (s = a_FirstSec; s < a_LastSec; s++)
	{
		// Print Info
		fprintf(LogOut, "CRL: Sector %u of %u\n", (uint32_t)s, (uint32_t)(a_LastSec - a_FirstSec));
		fflush(LogOut);
		
		// Render
		RenderTheView = 1;
		
		// Finish update (show to screen)
		// not rendering every single movement to the screen
		// speeds it up alot!
		if (RenderTheView)//if ((RenderRun & 511) == 0)
		{
			// Current status
			sprintf(Buf, "%i of %i", s, numsectors);
			M_WriteText(100, 100, Buf);
		
			// Show the avid watcher
			I_FinishUpdate();
		
			// Don't render anymore
			RenderTheView = 0;
		}
		
		// Get sector
		Sec = &sectors[s];
		
		// Get sector bounds
		sx = Sec->BBox[BOXLEFT];
		ex = Sec->BBox[BOXRIGHT];
		sy = Sec->BBox[BOXBOTTOM];
		ey = Sec->BBox[BOXTOP];
		
		// Round Coords
#if 1
		ddz = FIXED_TO_FLOAT(FRACSTEP);
		
		ddd = FIXED_TO_FLOAT(sx);
		ddd = trunc(((double)((int)(ddd / ddz))) * ddz);
		sx = FLOAT_TO_FIXED(ddd);
		
		ddd = FIXED_TO_FLOAT(sy);
		ddd = trunc(((double)((int)(ddd / ddz))) * ddz);
		sy = FLOAT_TO_FIXED(ddd);
		
		ddd = FIXED_TO_FLOAT(ex);
		ddd = trunc(((double)((int)(ddd / ddz))) * ddz);
		ex = FLOAT_TO_FIXED(ddd);
		
		ddd = FIXED_TO_FLOAT(ey);
		ddd = trunc(((double)((int)(ddd / ddz))) * ddz);
		ey = FLOAT_TO_FIXED(ddd);
#else
		sx = FixedMul(((FixedDiv(sx, FRACSTEP) >> FRACBITS) << FRACBITS), FRACSTEP);
		sy = FixedMul(((FixedDiv(sy, FRACSTEP) >> FRACBITS) << FRACBITS), FRACSTEP);
		ex = FixedMul(((FixedDiv(ex, FRACSTEP) >> FRACBITS) << FRACBITS), FRACSTEP);
		ey = FixedMul(((FixedDiv(ey, FRACSTEP) >> FRACBITS) << FRACBITS), FRACSTEP);
#endif
		
		// Correct?
		if (ex < sx)
		{
			z = ex;
			ex = sx;
			sx = z;
		}
		
		if (ey < sy)
		{
			z = ey;
			ey = sy;
			sy = z;
		}
		
		// Go through each coordinate (will take a very long time)
		for (x = sx + (FRACSTEP * a_BruteStepOffset); x < ex; x += (FRACSTEP * a_BruteStep))
			for (y = sy; y < ey; y += (FRACSTEP))
			{
				// Check if it isn't in the global bounding table
				if (x < bsx || x > bex || y < bsy || y > bey)
					continue;
					
				// Determine spot in image
#if 1
				ddz = FIXED_TO_FLOAT(FRACSTEP);
				
				ddd = ((FIXED_TO_FLOAT(x) - FIXED_TO_FLOAT(tMx)) / ddz);
				SpotX = trunc(ddd / a_BruteStep);
				
				ddd = ((FIXED_TO_FLOAT(y) - FIXED_TO_FLOAT(tMy)) / ddz);
				SpotY = trunc(ddd);
#else
				SpotX = (FixedDiv(x - tMx, FRACSTEP) >> FRACBITS) / a_BruteStep;
				SpotY = (FixedDiv(y - tMy, FRACSTEP) >> FRACBITS);
#endif
				// Sanity check
				if (SpotX < 0 || SpotY < 0 || SpotX >= ImW || SpotY >= ImH)
					continue;
				
				// Get current spot
				ThisHeat = &ImHeats[(ImW * SpotY) + SpotX];
				
				// See if subsector is valid
				SubS = R_IsPointInSubsector(x, y);
				RealSubS = R_PointInSubsector(x, y);
				
				// Nope!
				if (!SubS)
					continue;
				
				// Subsector is in another sector?
					// Don't want to render the same areas multiple times!
				//if (SubS->sector != Sec)
				//	continue;
			
				// Already checked?
				if (ThisHeat->NoWhere)
					continue;
				
				// Set as here
				ThisHeat->NoWhere = 1;
				ThisHeat->XPos = x >> FRACBITS;
				ThisHeat->YPos = y >> FRACBITS;
				
				// Put player subsector here (for heat map)
				players[0].mo->subsector = RealSubS;
				
				// Look around
				for (Ang = 0; Ang < 8; Ang++)
				{
					// Convert to angle_t
					RealAng = ANG45 * Ang;
					
					// Move player here
					players[0].mo->x = x;
					players[0].mo->y = y;
					players[0].mo->angle = RealAng;
					
					// 3 stage z position
					for (zzPos = 0; zzPos < 3; zzPos++, RenderRun++)
					{
						// Based on z position (floor, mid, ceil)
						if (zzPos == 0)
							players[0].mo->z = Sec->floorheight;
						else if (zzPos == 1)
							players[0].mo->z = Sec->floorheight + ((Sec->ceilingheight - Sec->floorheight) >> 1);
						else
							players[0].mo->z = Sec->ceilingheight;
						
						// Render
						PosChecks++;
						P_CalcHeight(&players[0]);
						R_RenderPlayerView(&players[0]);
						
						// Determine counts exceeds
							// Visplanes
						Temp = (int)(lastvisplane - visplanes);
						SawPlanes += Temp;
						if (Temp > ThisHeat->VisPlanes)
							ThisHeat->VisPlanes = Temp;
							// DrawSegs
						Temp = ActualDrawSegs;
						SawSegs += Temp;
						if (Temp > ThisHeat->DrawSegs)
							ThisHeat->DrawSegs = Temp;
							// Openings
						Temp = (lastopening - openings);
						SawOpenings += Temp;
						if (Temp > ThisHeat->Openings)
							ThisHeat->Openings = Temp;
							// Sprites
						Temp = ((vissprite_p - vissprites) + 1);
						SawSprites += Temp;
						if (Temp > ThisHeat->Sprites)
							ThisHeat->Sprites = Temp;
						
						// Automap activated? -- visplane overflow max
						if (automapactive)
						{
							AM_Stop();
							automapactive = 0;
						}
					}
				}
			}
	}
	
	/* Print Stats */
	fprintf(LogOut, "CRL: Checked %u positions\n", PosChecks);
	fprintf(LogOut, "CRL: Saw %u Planes\n", SawPlanes);
	fprintf(LogOut, "CRL: Saw %u Segs\n", SawSegs);
	fprintf(LogOut, "CRL: Saw %u Openings\n", SawOpenings);
	fprintf(LogOut, "CRL: Saw %u Sprites\n", SawSprites);
	fflush(LogOut);
	
	/* Dump heats to file/image pair */
	// Create files
	RawData = fopen((a_DATFile ? a_DATFile : "crlbf_rw.dat"), "wb");
	NicePic = fopen((a_PPMFile ? a_PPMFile : "crlbf_im.ppm"), "wb");
	
	// Dump it
	if (RawData || NicePic)
	{
		// Create headers for both
			// Raw Tables
		if (RawData)
		{
			u32 = ImW;
			fwrite(&u32, 4, 1, RawData);
			u32 = ImH;
			fwrite(&u32, 4, 1, RawData);
			u32 = tMx;
			fwrite(&u32, 4, 1, RawData);
			u32 = tMy;
			fwrite(&u32, 4, 1, RawData);
			u32 = FRACSTEP;
			fwrite(&u32, 4, 1, RawData);
		}
			// PPM
		fprintf(NicePic, "P6\n%i %i\n255\n", ImW, ImH * 4);
		
		// Dump raw data first
		if (RawData)
			for (i = 0; i < ImW * ImH; i++)
			{
				sixx = ImHeats[i].XPos;
				sixy = ImHeats[i].YPos;
				fwrite(&sixx, 2, 1, RawData);
				fwrite(&sixy, 2, 1, RawData);
	
				u8r = ImHeats[i].NoWhere;
				fwrite(&u8r, 1, 1, RawData);
			
				if (ImHeats[i].NoWhere)
				{
					u32 = ImHeats[i].VisPlanes;
					fwrite(&u32, 4, 1, RawData);
			
					u32 = ImHeats[i].DrawSegs;
					fwrite(&u32, 4, 1, RawData);
			
					u32 = ImHeats[i].Openings;
					fwrite(&u32, 4, 1, RawData);
			
					u32 = ImHeats[i].Sprites;
					fwrite(&u32, 4, 1, RawData);
				}
			}
		
		// Dump PPM now
		if (NicePic)
			for (j = 0; j < 4; j++)
				for (i = 0; i < ImW * ImH; i++)
				{
					// Dump which?
						// Nowhere on the map?
					if (!ImHeats[i].NoWhere)
					{
						u8r = 127;
						u8g = 127;
						if (j == 3)
							u8b = 0;
						else
							u8b = 127;
					}
					
						// Something is here
					else
					{
							// Visplanes
						if (j == 0)
						{
							u32 = ImHeats[i].VisPlanes;
							u32 = ((((double)u32) / (128.0 / 2)) * 255.0);
						
							// Really high (exceeds!)
							if (u32 > 511)
							{
								u8r = 255;
								u8g = 255;
								u8b = u32 - 511;
							}
						
							// Higher than half
							else if (u32 > 255)
							{
								u8r = u32 - 255;
								u8g = 255;
								u8b = 0;
							}
								// Less than zero?
							else if (u32 < 0)
							{
								u8r = 0;
								u8g = 0;
								u8b = 0;
							}
								// Normal
							else
							{
								u8r = 0;
								u8g = u32;
								u8b = 0;
							}
					
						}
							// DrawSegs
						else if (j == 1)
						{
							u32 = ImHeats[i].DrawSegs;
							u32 = ((((double)u32) / ((double)(NORMMAXDRAWSEGS / 2))) * 255.0);
						
							// Really high (exceeds!)
							if (u32 > 511)
							{
								u8r = 255;
								u8g = 255;
								u8b = u32 - 511;
							}
						
							// Higher than half
							else if (u32 > 255)
							{
								u8r = 255;
								u8g = u32 - 255;
								u8b = 0;
							}
								// Less than zero?
							else if (u32 < 0)
							{
								u8r = 0;
								u8g = 0;
								u8b = 0;
							}
								// Normal
							else
							{
								u8r = u32;
								u8g = 0;
								u8b = 0;
							}
						}
							// Openings
						else if (j == 2)
						{
							u32 = ImHeats[i].Openings;
							u32 = ((((double)u32) / ((double)(MAXOPENINGS / 2))) * 255.0);
					
							// Really high (exceeds!)
							if (u32 > 511)
							{
								u8r = 255;
								u8g = u32 - 511;
								u8b = 255;
							}
						
							// Higher than half
							else if (u32 > 255)
							{
								u8r = u32 - 255;
								u8g = 0;
								u8b = 255;
							}
								// Less than zero?
							else if (u32 < 0)
							{
								u8r = 0;
								u8g = 0;
								u8b = 0;
							}
								// Normal
							else
							{
								u8r = 0;
								u8g = 0;
								u8b = u32;
							}
						}
							// Sprites
						else if (j == 3)
						{
							u32 = ImHeats[i].Sprites;
							u32 = ((((double)u32) / ((double)(MAXVISSPRITES / 2))) * 255.0);
					
							// Really high (exceeds!)
							if (u32 > 511)
							{
								u8r = 255;
								u8g = 255;
								u8b = 255;
							}
						
							// Higher than half
							else if (u32 > 255)
							{
								u8r = 255;
								u8g = u32 - 255;
								u8b = 255;
							}
								// Less than zero?
							else if (u32 < 0)
							{
								u8r = 0;
								u8g = 0;
								u8b = 0;
							}
								// Normal
							else
							{
								u8r = u32;
								u8g = u32;
								u8b = u32;
							}
						}
					}
			
					// Write
					fwrite(&u8r, 1, 1, NicePic);
					fwrite(&u8g, 1, 1, NicePic);
					fwrite(&u8b, 1, 1, NicePic);
				}
		
		// Close files
		if (RawData)
			fclose(RawData);
		if (NicePic)
			fclose(NicePic);
	}
	
	fprintf(LogOut, "CRL: Wrote Files!\n");
	fflush(LogOut);
	
	// Close log
	if (LogOut != stdout)
		fclose(LogOut);
	
	// Free it
	free(ImHeats);
	
	/* Set player position */
	players[0].mo->x = OldX;
	players[0].mo->y = OldY;
	players[0].mo->z = OldZ;
	P_SetThingPosition(players[0].mo);
}


