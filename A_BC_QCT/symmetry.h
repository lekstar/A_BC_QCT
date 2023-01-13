#ifndef SYMMETRY_H
#define SYMMETRY_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795028841971694
#endif

#define	DIMENSION 3
#define MAXPARAM  7

    typedef struct {
        int     type;
        double  x[DIMENSION];
    } ATOM1;

    typedef struct _SYMMETRY_ELEMENT_ {
        void    (*transform_atom)(struct _SYMMETRY_ELEMENT_* el, ATOM1* from, ATOM1* to);
        int* transform;     /*   Correspondence table for the transformation         */
        int     order;         /*   Applying transformation this many times is identity */
        int     nparam;        /*   4 for inversion and planes, 7 for axes              */
        double  maxdev;        /*   Larges error associated with the element            */
        double  distance;
        double  normal[DIMENSION];
        double  direction[DIMENSION];
        int	rank;		/*   Number of unique atoms for this element		  */
    } SYMMETRY_ELEMENT;

    typedef struct {
        char* group_name;        /* Canonical group name                              */
        char* symmetry_code;     /* Group symmetry code                               */
        int     (*check)(void);  /* Additional verification routine, not used         */
    } POINT_GROUP;

    extern	double                 ToleranceSame;
    extern	double                 TolerancePrimary;
    extern	double                 ToleranceFinal;
    extern	double                 MaxOptStep;
    extern	double                 MinOptStep;
    extern	double                 GradientStep;
    extern	double                 OptChangeThreshold;
    extern	double                 CenterOfSomething[DIMENSION];
    extern	double* DistanceFromCenter;
    extern	int                    verbose;
    extern	int                    MaxOptCycles;
    extern	int                    OptChangeHits;
    extern	int                    MaxAxisOrder;
    extern	int                    AtomsCount;
    extern	ATOM1* Atoms;
    extern	int                    PlanesCount;
    extern	SYMMETRY_ELEMENT** Planes;
    extern	SYMMETRY_ELEMENT* MolecularPlane;
    extern	int                    InversionCentersCount;
    extern	SYMMETRY_ELEMENT** InversionCenters;
    extern	int                    NormalAxesCount;
    extern	SYMMETRY_ELEMENT** NormalAxes;
    extern	int                    ImproperAxesCount;
    extern	SYMMETRY_ELEMENT** ImproperAxes;
    extern	int* NormalAxesCounts;
    extern	int* ImproperAxesCounts;
    extern	int                    BadOptimization;
    extern	char* SymmetryCode;
    extern	int                    OrdUniq;
    extern	int                    OrdAxis;
    extern	int                    SymmPlane;
    extern	int                    SymmAxis;
    extern	int                    SymmImproper;
    extern	int                    SymmCenter;
    extern	int                    NonAbel;
    /*
     *    Statistics
     */
    extern	long                   StatTotal;
    extern	long                   StatEarly;
    extern	long                   StatPairs;
    extern	long                   StatDups;
    extern	long                   StatOrder;
    extern	long                   StatOpt;
    extern	long                   StatAccept;

#ifdef __cplusplus
    extern "C" {
#endif

    int point_group_calc(int num_atom, double* charge, double* coord, char* PG, double fin_tol);

#ifdef __cplusplus
    }
#endif

#endif
