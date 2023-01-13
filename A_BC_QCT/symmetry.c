/*
 *  Brute force symmetry analyzer.
 *  This is actually C++ program, masquerading as a C one!
 *
 *  (C) 1996 S. Pachkovsky, ps@oci.unizh.ch
 *  This code comes with NO WARRANTY, neither explicit nor implied.
 *  Non-profit use of the code is welcome, provided that the
 *  copyright above is retained in all derived works and the
 *  origin of the code is clearly indicated. Any for-profit use
 *  requires prior written permission of the author.
 *
 *	(C) 2002-2004 Georgy Salnikov, sge@nmr.nioch.nsc.ru, 20040324 SGE
 *
 *	The code enhanced to provide the following new features:
 *
 *	Reorienting the symmetry elements and the molecule itself
 *	so that the main (abelian or not) symmetry elements
 *	be aligned along the master frame coordinate axes
 *	Reoriented symmetry elements and coordinates of atoms are printed out
 *
 *	Symmetrizing the atoms according to all or a single symmetry element
 *	in order to eliminate possible minor symmetry distortions
 *	Corrected coordinates of all atoms and of symmetry unique atoms
 *	(according to all symmetry elements) are printed out
 *
 *	Abelian symmetry elements, symmetry generators and symmetry subgroup
 *	detected and printed out
 *
 *	Symmetrizing the atoms according to abelian symmetry elements
 *	Corrected coordinates of all atoms and of symmetry unique atoms
 *	(according to abelian symmetry elements only) are printed out
 *
 * $Log: symmetry.c,v $
 * Revision 1.15  2000/01/25  16:47:17  patchkov
 * *** empty log message ***
 *
 * Revision 1.14  2000/01/25  16:39:08  patchkov
 * *** empty log message ***
 *
 * Revision 1.13  1996/05/24  12:32:08  ps
 * Revision 1.12  1996/05/23  16:10:47  ps
 * First reasonably stable version.
 *
 */

#include "symmetry.h"

/*
 *  All specific structures should have corresponding elements in the
 *  same position generic structure does.
 *
 *  Planes are characterized by the surface normal direction
 *  (taken in the direction *from* the coordinate origin)
 *  and distance from the coordinate origin to the plane
 *  in the direction of the surface normal.
 *
 *  Inversion is characterized by location of the inversion center.
 *
 *  Rotation is characterized by a vector (distance+direction) from the origin
 *  to the rotation axis, axis direction and rotation order. Rotations
 *  are in the clockwise direction looking opposite to the direction
 *  of the axis. Note that this definition of the rotation axis
 *  is *not* unique, since an arbitrary multiple of the axis direction
 *  can be added to the position vector without changing actual operation.
 *
 *  Mirror rotation is defined by the same parameters as normal rotation,
 *  but the origin is now unambiguous since it defines the position of the
 *  plane associated with the axis.
 *
 */



double                 ToleranceSame         = 1e-3 ;
double                 TolerancePrimary      = 5e-2 ;
double                 ToleranceFinal        = 1e-4 ;
double                 MaxOptStep            = 5e-1 ;
double                 MinOptStep            = 1e-7 ;
double                 GradientStep          = 1e-7 ;
double                 OptChangeThreshold    = 1e-10 ;
double                 CenterOfSomething[ DIMENSION ] ;
double *               DistanceFromCenter    = NULL ;
int                    verbose               = -3 ;
int                    MaxOptCycles          = 1000 ;
int                    OptChangeHits         = 5 ;
int                    MaxAxisOrder          = 20 ;
int                    AtomsCount            = 0 ;
ATOM1 *                 Atoms                 = NULL ;
int                    PlanesCount           = 0 ;
SYMMETRY_ELEMENT **    Planes                = NULL ;
SYMMETRY_ELEMENT *     MolecularPlane        = NULL ;
int                    InversionCentersCount = 0 ;
SYMMETRY_ELEMENT **    InversionCenters      = NULL ;
int                    NormalAxesCount       = 0 ;
SYMMETRY_ELEMENT **    NormalAxes            = NULL ;
int                    ImproperAxesCount     = 0 ;
SYMMETRY_ELEMENT **    ImproperAxes          = NULL ;
int *                  NormalAxesCounts      = NULL ;
int *                  ImproperAxesCounts    = NULL ;
int                    BadOptimization       = 0 ;
char *                 SymmetryCode          = "" ;
int                    OrdUniq               = 1 ;
int                    OrdAxis               = 1 ;
int                    SymmPlane             = 0 ;
int                    SymmAxis              = 0 ;
int                    SymmImproper          = 0 ;
int                    SymmCenter            = 0 ;
int                    NonAbel               = 0 ;
/*
 *    Statistics
 */
long                   StatTotal             = 0 ;
long                   StatEarly             = 0 ;
long                   StatPairs             = 0 ;
long                   StatDups              = 0 ;
long                   StatOrder             = 0 ;
long                   StatOpt               = 0 ;
long                   StatAccept            = 0 ;

/*
 *    Point groups I know about
 */
int true(void){ return 1 ; }
POINT_GROUP            PointGroups[]         = {
    {  "C1",    "",                                                          true  },
    {  "Cs",    "(sigma) ",                                                  true  },
    {  "Ci",    "(i) ",                                                      true  },
    {  "C2",    "(C2) ",                                                     true  },
    {  "C3",    "(C3) ",                                                     true  },
    {  "C4",    "(C4) (C2) ",                                                true  },
    {  "C5",    "(C5) ",                                                     true  },
    {  "C6",    "(C6) (C3) (C2) ",                                           true  },
    {  "C7",    "(C7) ",                                                     true  },
    {  "C8",    "(C8) (C4) (C2) ",                                           true  },
    {  "D2",    "3*(C2) ",                                                   true  },
    {  "D3",    "(C3) 3*(C2) ",                                              true  },
    {  "D4",    "(C4) 5*(C2) ",                                              true  },
    {  "D5",    "(C5) 5*(C2) ",                                              true  },
    {  "D6",    "(C6) (C3) 7*(C2) ",                                         true  },
    {  "D7",    "(C7) 7*(C2) ",                                              true  },
    {  "D8",    "(C8) (C4) 9*(C2) ",                                         true  },
    {  "C2v",   "(C2) 2*(sigma) ",                                           true  },
    {  "C3v",   "(C3) 3*(sigma) ",                                           true  },
    {  "C4v",   "(C4) (C2) 4*(sigma) ",                                      true  },
    {  "C5v",   "(C5) 5*(sigma) ",                                           true  },
    {  "C6v",   "(C6) (C3) (C2) 6*(sigma) ",                                 true  },
    {  "C7v",   "(C7) 7*(sigma) ",                                           true  },
    {  "C8v",   "(C8) (C4) (C2) 8*(sigma) ",                                 true  },
    {  "C2h",   "(i) (C2) (sigma) ",                                         true  },
    {  "C3h",   "(C3) (S3) (sigma) ",                                        true  },
    {  "C4h",   "(i) (C4) (C2) (S4) (sigma) ",                               true  },
    {  "C5h",   "(C5) (S5) (sigma) ",                                        true  },
    {  "C6h",   "(i) (C6) (C3) (C2) (S6) (S3) (sigma) ",                     true  },
    {  "C7h",   "(C7) (S7) (sigma) ",                                        true  },
    {  "C8h",   "(i) (C8) (C4) (C2) (S8) (S4) (sigma) ",                     true  },
    {  "D2h",   "(i) 3*(C2) 3*(sigma) ",                                     true  },
    {  "D3h",   "(C3) 3*(C2) (S3) 4*(sigma) ",                               true  },
    {  "D4h",   "(i) (C4) 5*(C2) (S4) 5*(sigma) ",                           true  },
    {  "D5h",   "(C5) 5*(C2) (S5) 6*(sigma) ",                               true  },
    {  "D6h",   "(i) (C6) (C3) 7*(C2) (S6) (S3) 7*(sigma) ",                 true  },
    {  "D7h",   "(C7) 7*(C2) (S7) 8*(sigma) ",                               true  },
    {  "D8h",   "(i) (C8) (C4) 9*(C2) (S8) (S4) 9*(sigma) ",                 true  },
    {  "D2d",   "3*(C2) (S4) 2*(sigma) ",                                    true  },
    {  "D3d",   "(i) (C3) 3*(C2) (S6) 3*(sigma) ",                           true  },
    {  "D4d",   "(C4) 5*(C2) (S8) 4*(sigma) ",                               true  },
    {  "D5d",   "(i) (C5) 5*(C2) (S10) 5*(sigma) ",                          true  },
    {  "D6d",   "(C6) (C3) 7*(C2) (S12) (S4) 6*(sigma) ",                    true  },
    {  "D7d",   "(i) (C7) 7*(C2) (S14) 7*(sigma) ",                          true  },
    {  "D8d",   "(C8) (C4) 9*(C2) (S16) 8*(sigma) ",                         true  },
    {  "S4",    "(C2) (S4) ",                                                true  },
    {  "S6",    "(i) (C3) (S6) ",                                            true  },
    {  "S8",    "(C4) (C2) (S8) ",                                           true  },
    {  "T",     "4*(C3) 3*(C2) ",                                            true  },
    {  "Th",    "(i) 4*(C3) 3*(C2) 4*(S6) 3*(sigma) ",                       true  },
    {  "Td",    "4*(C3) 3*(C2) 3*(S4) 6*(sigma) ",                           true  },
    {  "O",     "3*(C4) 4*(C3) 9*(C2) ",                                     true  },
    {  "Oh",    "(i) 3*(C4) 4*(C3) 9*(C2) 4*(S6) 3*(S4) 9*(sigma) ",         true  },
    {  "Cinfv", "(Cinf) (sigma) ",                                           true  },
    {  "Dinfh", "(i) (Cinf) (C2) 2*(sigma) ",                                true  },
    {  "I",     "6*(C5) 10*(C3) 15*(C2) ",                                   true  },
    {  "Ih",    "(i) 6*(C5) 10*(C3) 15*(C2) 6*(S10) 10*(S6) 15*(sigma) ",    true  },
    {  "Kh",    "(i) (Cinf) (sigma) ",                                       true  },
    } ;
#define PointGroupsCount (sizeof(PointGroups)/sizeof(POINT_GROUP))
char *                 PointGroupRejectionReason = NULL ;

/*
 *   Generic functions
 */

double
pow2( double x )
{
return x * x ;
}

int
establish_pairs( SYMMETRY_ELEMENT *elem )
{
        int               i, j, k, best_j ;
        char *            atom_used = calloc( AtomsCount, 1 ) ;
        double            distance, best_distance ;
        ATOM1              symmetric ;

if( atom_used == NULL ){
    fprintf( stderr, "Out of memory for tagging array in establish_pairs()\n" ) ;
    exit( EXIT_FAILURE ) ;
    }
for( i = 0 ; i < AtomsCount ; i++ ){
    if( elem->transform[i] >= AtomsCount ){ /* No symmetric atom yet          */
        if( verbose > 2 ) printf( "        looking for a pair for %d\n", i ) ;
        elem->transform_atom( elem, Atoms+i, &symmetric ) ;
        if( verbose > 2 ) printf( "        new coordinates are: (%g,%g,%g)\n", 
                              symmetric.x[0], symmetric.x[1], symmetric.x[2] ) ;
        best_j        = i ;
        best_distance = 2*TolerancePrimary ;/* Performance value we'll reject */
        for( j = 0 ; j < AtomsCount ; j++ ){
            if( Atoms[j].type != symmetric.type || atom_used[j] )
                continue ;
            for( k = 0, distance = 0 ; k < DIMENSION ; k++ ){
                distance += pow2( symmetric.x[k] - Atoms[j].x[k] ) ;
                }
            distance = sqrt( distance ) ;
            if( verbose > 2 ) printf( "        distance to %d is %g\n", j, distance ) ;
            if( distance < best_distance ){
                best_j        = j ;
                best_distance = distance ;
                }
            }
        if( best_distance > TolerancePrimary ){ /* Too bad, there is no symmetric atom */
            if( verbose > 0 ) 
                printf( "        no pair for atom %d - best was %d with err = %g\n", i, best_j, best_distance ) ;
            free( atom_used ) ;
            return -1 ;
            }
        elem->transform[i] = best_j ;
        atom_used[best_j]  = 1 ;
        if( verbose > 1 ) printf( "        atom %d transforms to the atom %d, err = %g\n", i, best_j, best_distance ) ;
        }
    }
free( atom_used ) ;
return 0 ;
}

int
check_transform_order( SYMMETRY_ELEMENT *elem )
{
        int             i, j, k ;
        void            rotate_reflect_atom( SYMMETRY_ELEMENT *, ATOM1 *, ATOM1 *) ;

for( i = 0 ; i < AtomsCount ; i++ ){
    if( elem->transform[i] == i )   /* Identity transform is Ok for any order */
        continue ;
    if( elem->transform_atom == rotate_reflect_atom ){
        j = elem->transform[i] ;
        if( elem->transform[j] == i )
            continue ; /* Second-order transform is Ok for improper axis */
        }
    for( j = elem->order - 1, k = elem->transform[i] ; j > 0 ; j--, k = elem->transform[k] ){
        if( k == i ){
            if( verbose > 0 ) printf( "        transform looped %d steps too early from atom %d\n", j, i ) ;
            return -1 ;
            }
        }
    if( k != i && elem->transform_atom == rotate_reflect_atom ){
        /* For improper axes, the complete loop may also take twice the order */
        for( j = elem->order ; j > 0 ; j--, k = elem->transform[k] ){
            if( k == i ){
                if( verbose > 0 ) printf( "        (improper) transform looped %d steps too early from atom %d\n", j, i ) ;
                return -1 ;
                }
            }
        }
    if( k != i ){
        if( verbose > 0 ) printf( "        transform failed to loop after %d steps from atom %d\n", elem->order, i ) ;
        return -1 ;
        }
    }
return 0 ;
}

int
same_transform( SYMMETRY_ELEMENT *a, SYMMETRY_ELEMENT *b )
{
        int               i, j ;
        int               code ;

if( ( a->order != b->order ) || ( a->nparam != b->nparam ) || ( a->transform_atom != b->transform_atom ) )
    return 0 ;
for( i = 0, code = 1 ; i < AtomsCount ; i++ ){
    if( a->transform[i] != b->transform[i] ){
        code = 0 ;
        break ;
        }
    }
if( code == 0 && a->order > 2 ){  /* b can also be a reverse transformation for a */
    for( i = 0 ; i < AtomsCount ; i++ ){
        j = a->transform[i] ;
        if( b->transform[j] != i )
            return 0 ;
        }
    return 1 ;
    }
return code ;
}

SYMMETRY_ELEMENT *
alloc_symmetry_element( void )
{
        SYMMETRY_ELEMENT * elem = calloc( 1, sizeof( SYMMETRY_ELEMENT ) ) ;
        int                i ;

if( elem == NULL ){
    fprintf( stderr, "Out of memory allocating symmetry element\n" ) ;
    exit( EXIT_FAILURE ) ;
    }
elem->transform = calloc( AtomsCount, sizeof( int ) ) ;
if( elem->transform == NULL ){
    fprintf( stderr, "Out of memory allocating transform table for symmetry element\n" ) ;
    exit( EXIT_FAILURE ) ;
    }
for( i = 0 ; i < AtomsCount ; i++ ){
    elem->transform[i] = AtomsCount + 1 ; /* An impossible value */
    }
return elem ;
}

void
destroy_symmetry_element( SYMMETRY_ELEMENT *elem )
{
if( elem != NULL ){
    if( elem->transform != NULL )
        free( elem->transform ) ;
    free( elem ) ;
    }
}

int
check_transform_quality( SYMMETRY_ELEMENT *elem )
{
        int               i, j, k ;
        ATOM1              symmetric ;
        double            r, max_r ;

for( i = 0, max_r = 0 ; i < AtomsCount ; i++ ){
    j = elem->transform[i] ;
    elem->transform_atom( elem, Atoms + i, &symmetric ) ;
    for( k = 0, r = 0 ; k < DIMENSION ; k++ ){
        r += pow2( symmetric.x[k] - Atoms[j].x[k] ) ;
        }
    r = sqrt( r ) ;
    if( r > ToleranceFinal ){
        if( verbose > 0 ) printf( "        distance to symmetric atom (%g) is too big for %d\n", r, i ) ;
        return -1 ;
        }
    if( r > max_r ) max_r = r ;
    }
elem->maxdev = max_r ;
return 0 ;
}

double
eval_optimization_target_function( SYMMETRY_ELEMENT *elem, int *finish )
{
        int               i, j, k ;
        ATOM1              symmetric ;
        double            target, r, maxr ;

if( elem->nparam >= 4 ){
    for( k = 0, r = 0 ; k < DIMENSION ; k++ ){
        r += elem->normal[k]*elem->normal[k] ;
        }
    r = sqrt( r ) ;
    if( r < ToleranceSame ){
        fprintf( stderr, "Normal collapced!\n" ) ;
        exit( EXIT_FAILURE ) ;
        }
    for( k = 0 ; k < DIMENSION ; k++ ){
        elem->normal[k] /= r ;
        }
    if( elem->distance < 0 ){
        elem->distance = -elem->distance ;
        for( k = 0 ; k < DIMENSION ; k++ ){
            elem->normal[k] = -elem->normal[k] ;
            }
        }
    }
if( elem->nparam >= 7 ){
    for( k = 0, r = 0 ; k < DIMENSION ; k++ ){
        r += elem->direction[k]*elem->direction[k] ;
        }
    r = sqrt( r ) ;
    if( r < ToleranceSame ){
        fprintf( stderr, "Direction collapced!\n" ) ;
        exit( EXIT_FAILURE ) ;
        }
    for( k = 0 ; k < DIMENSION ; k++ ){
        elem->direction[k] /= r ;
        }
    }
for( i = 0, target = maxr = 0 ; i < AtomsCount ; i++ ){
    elem->transform_atom( elem, Atoms + i, &symmetric ) ;
    j = elem->transform[i] ;
    for( k = 0, r = 0 ; k < DIMENSION ; k++ ){
        r += pow2( Atoms[j].x[k] - symmetric.x[k] ) ;
        }
    if( r > maxr ) maxr = r ;
    target += r ;
    }
if( finish != NULL ){
    *finish = 0 ;
    if( sqrt( maxr ) < ToleranceFinal )
        *finish = 1 ;
    }
return target ;
}

void
get_params( SYMMETRY_ELEMENT *elem, double values[] )
{
memcpy( values, &elem->distance, elem->nparam * sizeof( double ) ) ;
}

void
set_params( SYMMETRY_ELEMENT *elem, double values[] )
{
memcpy( &elem->distance, values, elem->nparam * sizeof( double ) ) ;
}

void
optimize_transformation_params( SYMMETRY_ELEMENT *elem )
{
        double            values[ MAXPARAM ] ;
        double            grad  [ MAXPARAM ] ;
        double            force [ MAXPARAM ] ;
        double            step  [ MAXPARAM ] ;
        double            f, fold, fnew, fnew2, fdn, fup, snorm ;
        double            a, b, x ;
        int               vars  = elem->nparam ;
        int               cycle = 0 ;
        int               i, finish ;
        int               hits = 0 ;

if( vars > MAXPARAM ){
    fprintf( stderr, "Catastrophe in optimize_transformation_params()!\n" ) ;
    exit( EXIT_FAILURE ) ;
    }
f = 0 ;
do {
    fold = f ;
    f    = eval_optimization_target_function( elem, &finish ) ;
    /* Evaluate function, gradient and diagonal force constants */
    if( verbose > 1 ) printf( "            function value = %g\n", f ) ;
    if( finish ){
        if( verbose > 1 ) printf( "        function value is small enough\n" ) ;
        break ;
        }
    if( cycle > 0 ){
        if( fabs( f-fold ) > OptChangeThreshold )
             hits = 0 ;
        else hits++ ;
        if( hits >= OptChangeHits ){
            if( verbose > 1 ) printf( "        no progress is made, stop optimization\n" ) ;
            break ;
            }
        }
    get_params( elem, values ) ;
    for( i = 0 ; i < vars ; i++ ){
        values[i] -= GradientStep ;
        set_params( elem, values ) ;
        fdn        = eval_optimization_target_function( elem, NULL ) ;
        values[i] += 2*GradientStep ;
        set_params( elem, values ) ;
        fup        = eval_optimization_target_function( elem, NULL ) ;
        values[i] -= GradientStep ;
        grad[i]    = ( fup - fdn ) / ( 2 * GradientStep ) ;
        force[i]   = ( fup + fdn - 2*f ) / ( GradientStep * GradientStep ) ;
        if( verbose > 1 ) printf( "        i = %d, grad = %12.6e, force = %12.6e\n", i, grad[i], force[i] ) ;
        }
    /* Do a quasy-Newton step */
    for( i = 0, snorm = 0 ; i < vars ; i++ ){
        if( force[i] <  0   ) force[i] = -force[i] ;
        if( force[i] < 1e-3 ) force[i] = 1e-3 ;
        if( force[i] > 1e3  ) force[i] = 1e3 ;
        step[i] = - grad[i]/force[i] ;
        snorm += step[i] * step[i] ;
        }
    snorm = sqrt( snorm ) ;
    if( snorm > MaxOptStep ){ /* Renormalize step */
        for( i = 0 ; i < vars ; i++ )
            step[i] *= MaxOptStep/snorm ;
        snorm = MaxOptStep ;
        }
    do {
        for( i = 0 ; i < vars ; i++ ){
            values[i] += step[i] ;
            }
        set_params( elem, values ) ;
        fnew = eval_optimization_target_function( elem, NULL ) ;
        if( fnew < f )
            break ;
        for( i = 0 ; i < vars ; i++ ){
            values[i] -= step[i] ;
            step  [i] /= 2 ;
            }
        set_params( elem, values ) ;
        snorm /= 2 ;
        } while( snorm > MinOptStep ) ;
        if( (snorm > MinOptStep) && (snorm < MaxOptStep / 2) ){  /* try to do quadratic interpolation */
            for( i = 0 ; i < vars ; i++ )
                values[i] += step[i] ;
            set_params( elem, values ) ;
            fnew2 = eval_optimization_target_function( elem, NULL ) ;
            if( verbose > 1 ) printf( "        interpolation base points: %g, %g, %g\n", f, fnew, fnew2 ) ;
            for( i = 0 ; i < vars ; i++ )
                values[i] -= 2*step[i] ;
            a     = ( 4*f - fnew2 - 3*fnew ) / 2 ;
            b     = ( f + fnew2 - 2*fnew ) / 2 ;
            if( verbose > 1 ) printf( "        linear interpolation coefficients %g, %g\n", a, b ) ;
            if( b > 0 ){
                x = -a/(2*b) ;
                if( x > 0.2 && x < 1.8 ){
                    if( verbose > 1 ) printf( "        interpolated: %g\n", x ) ;
                    for( i = 0 ; i < vars ; i++ )
                        values[i] += x*step[i] ;
                    }
                else b = 0 ;
                }
            if( b <= 0 ){
                if( fnew2 < fnew ){
                    for( i = 0 ; i < vars ; i++ )
                        values[i] += 2*step[i] ;
                    }
                else {
                    for( i = 0 ; i < vars ; i++ )
                        values[i] += step[i] ;
                    }
                }
            set_params( elem, values ) ;
            }
    } while( snorm > MinOptStep && ++cycle < MaxOptCycles ) ;
f = eval_optimization_target_function( elem, NULL ) ;
if( cycle >= MaxOptCycles ) BadOptimization = 1 ;
if( verbose > 0 ) {
    if( cycle >= MaxOptCycles )
        printf( "        maximum number of optimization cycles made\n" ) ;
        printf( "        optimization completed after %d cycles with f = %g\n", cycle, f ) ;
    }
}

int
refine_symmetry_element( SYMMETRY_ELEMENT *elem, int build_table )
{
        int               i ;


if( build_table && (establish_pairs( elem ) < 0) ){
    StatPairs++ ;
    if( verbose > 0 ) printf( "        no transformation correspondence table can be constructed\n" ) ;
    return -1 ;
    }
for( i = 0 ; i < PlanesCount ; i++ ){
    if( same_transform( Planes[i], elem ) ){
        StatDups++ ;
        if( verbose > 0 ) printf( "        transformation is identical to plane %d\n", i ) ;
        return -1 ;
        }
    }
for( i = 0 ; i < InversionCentersCount ; i++ ){
    if( same_transform( InversionCenters[i], elem ) ){
        StatDups++ ;
        if( verbose > 0 ) printf( "        transformation is identical to inversion center %d\n", i ) ;
        return -1 ;
        }
    }
for( i = 0 ; i < NormalAxesCount ; i++ ){
    if( same_transform( NormalAxes[i], elem ) ){
        StatDups++ ;
        if( verbose > 0 ) printf( "        transformation is identical to normal axis %d\n", i ) ;
        return -1 ;
        }
    }
for( i = 0 ; i < ImproperAxesCount ; i++ ){
    if( same_transform( ImproperAxes[i], elem ) ){
        StatDups++ ;
        if( verbose > 0 ) printf( "        transformation is identical to improper axis %d\n", i ) ;
        return -1 ;
        }
    }
if( check_transform_order( elem ) < 0 ){
    StatOrder++ ;
    if( verbose > 0 ) printf( "        incorrect transformation order\n" ) ;
    return -1 ;
    }
optimize_transformation_params( elem ) ;
if( check_transform_quality( elem ) < 0 ){
    StatOpt++ ;
    if( verbose > 0 ) printf( "        refined transformation does not pass the numeric threshold\n" ) ;
    return -1 ;
    }
StatAccept++ ;
return 0 ;
}

/*
 *   Plane-specific functions
 */

void
mirror_atom( SYMMETRY_ELEMENT *plane, ATOM1 *from, ATOM1 *to )
{
        int                i ;
        double             r ;

for( i = 0, r = plane->distance ; i < DIMENSION ; i++ ){
    r -= from->x[i] * plane->normal[i] ;
    }
to->type = from->type ;
for( i = 0 ; i < DIMENSION ; i++ ){
    to->x[i] = from->x[i] + 2*r*plane->normal[i] ;
    }
}

SYMMETRY_ELEMENT *
init_mirror_plane( int i, int j )
{
        SYMMETRY_ELEMENT * plane = alloc_symmetry_element() ;
        double             dx[ DIMENSION ], midpoint[ DIMENSION ], rab, r ;
        int                k ;

if( verbose > 0 ) printf( "Trying mirror plane for atoms %d,%d\n", i, j ) ;
StatTotal++ ;
plane->transform_atom = mirror_atom ;
plane->order          = 2 ;
plane->nparam         = 4 ;
for( k = 0, rab = 0 ; k < DIMENSION ; k++ ){
    dx[k]       = Atoms[i].x[k] - Atoms[j].x[k] ;
    midpoint[k] = ( Atoms[i].x[k] + Atoms[j].x[k] ) / 2.0 ;
    rab        += dx[k]*dx[k] ;
    }
rab = sqrt(rab) ;
if( rab < ToleranceSame ){
    fprintf( stderr, "Atoms %d and %d coincide (r = %g)\n", i, j, rab ) ;
    exit( EXIT_FAILURE ) ;
    }
for( k = 0, r = 0 ; k < DIMENSION ; k++ ){
    plane->normal[k] = dx[k]/rab ;
    r += midpoint[k]*plane->normal[k] ;
    }
if( r < 0 ){  /* Reverce normal direction, distance is always positive! */
    r = -r ;
    for( k = 0 ; k < DIMENSION ; k++ ){
        plane->normal[k] = -plane->normal[k] ;
        }
    }
plane->distance = r ;
if( verbose > 0 ) printf( "    initial plane is at %g from the origin\n", r ) ;
if( refine_symmetry_element( plane, 1 ) < 0 ){
    if( verbose > 0 ) printf( "    refinement failed for the plane\n" ) ;
    destroy_symmetry_element( plane ) ;
    return NULL ;
    }
return plane ;
}

SYMMETRY_ELEMENT *
init_ultimate_plane( void )
{
        SYMMETRY_ELEMENT * plane = alloc_symmetry_element() ;
        double             d0[ DIMENSION ], d1[ DIMENSION ], d2[ DIMENSION ] ;
        double             p[ DIMENSION ] ;
        double             r, s0, s1, s2 ;
        double *           d ;
        int                i, j, k ;

if( verbose > 0 ) printf( "Trying whole-molecule mirror plane\n" ) ;
StatTotal++ ;
plane->transform_atom = mirror_atom ;
plane->order          = 1 ;
plane->nparam         = 4 ;
for( k = 0 ; k < DIMENSION ; k++ )
    d0[k] = d1[k] = d2[k] = 0 ;
d0[0] = 1 ; d1[1] = 1 ; d2[2] = 1 ;
for( i = 1 ; i < AtomsCount ; i++ ){
    for( j = 0 ; j < i ; j++ ){
        for( k = 0, r = 0 ; k < DIMENSION ; k++ ){
            p[k] = Atoms[i].x[k] - Atoms[j].x[k] ;
            r   += p[k]*p[k] ;
            }
        r = sqrt(r) ;
        for( k = 0, s0=s1=s2=0 ; k < DIMENSION ; k++ ){
            p[k] /= r ;
            s0   += p[k]*d0[k] ;
            s1   += p[k]*d1[k] ;
            s2   += p[k]*d2[k] ;
            }
        for( k = 0 ; k < DIMENSION ; k++ ){
            d0[k] -= s0*p[k] ;
            d1[k] -= s1*p[k] ;
            d2[k] -= s2*p[k] ;
            }
        }
    }
for( k = 0, s0=s1=s2=0 ; k < DIMENSION ; k++ ){
    s0 += d0[k] ;
    s1 += d1[k] ;
    s2 += d2[k] ;
    }
d = NULL ;
if( s0 >= s1 && s0 >= s2 ) d = d0 ;
if( s1 >= s0 && s1 >= s2 ) d = d1 ;
if( s2 >= s0 && s2 >= s1 ) d = d2 ;
if( d == NULL ){
    fprintf( stderr, "Catastrophe in init_ultimate_plane(): %g, %g and %g have no ordering!\n", s0, s1, s2 ) ;
    exit( EXIT_FAILURE ) ;
    }
for( k = 0, r = 0 ; k < DIMENSION ; k++ )
    r += d[k]*d[k] ;
r = sqrt(r) ;
if( r > 0 ){
    for( k = 0 ; k < DIMENSION ; k++ )
        plane->normal[k] = d[k]/r ;
    }
else {
    for( k = 1 ; k < DIMENSION ; k++ )
        plane->normal[k] = 0 ;
    plane->normal[0] = 1 ;
    }
for( k = 0, r = 0 ; k < DIMENSION ; k++ )
    r += CenterOfSomething[k]*plane->normal[k] ;
plane->distance = r ;
for( k = 0 ; k < AtomsCount ; k++ )
    plane->transform[k] = k ;
if( refine_symmetry_element( plane, 0 ) < 0 ){
    if( verbose > 0 ) printf( "    refinement failed for the plane\n" ) ;
    destroy_symmetry_element( plane ) ;
    return NULL ;
    }
return plane ;
}
/*
 *   Inversion-center specific functions
 */
void
invert_atom( SYMMETRY_ELEMENT *center, ATOM1 *from, ATOM1 *to )
{
        int                i ;

to->type = from->type ;
for( i = 0 ; i < DIMENSION ; i++ ){
    to->x[i] = 2*center->distance*center->normal[i] - from->x[i] ;
    }
}

SYMMETRY_ELEMENT *
init_inversion_center( void )
{
        SYMMETRY_ELEMENT * center = alloc_symmetry_element() ;
        int                k ;
        double             r ;

if( verbose > 0 ) printf( "Trying inversion center at the center of something\n" ) ;
StatTotal++ ;
center->transform_atom = invert_atom ;
center->order          = 2 ;
center->nparam         = 4 ;
for( k = 0, r = 0 ; k < DIMENSION ; k++ )
    r += CenterOfSomething[k]*CenterOfSomething[k] ;
r = sqrt(r) ;
if( r > 0 ){
    for( k = 0 ; k < DIMENSION ; k++ )
        center->normal[k] = CenterOfSomething[k]/r ;
    }
else {
    center->normal[0] = 1 ;
    for( k = 1 ; k < DIMENSION ; k++ )
        center->normal[k] = 0 ;
    }
center->distance = r ;
if( verbose > 0 ) printf( "    initial inversion center is at %g from the origin\n", r ) ;
if( refine_symmetry_element( center, 1 ) < 0 ){
    if( verbose > 0 ) printf( "    refinement failed for the inversion center\n" ) ;
    destroy_symmetry_element( center ) ;
    return NULL ;
    }
return center ;
}

/*
 *   Normal rotation axis-specific routines.
 */
void
rotate_atom( SYMMETRY_ELEMENT *axis, ATOM1 *from, ATOM1 *to )
{
        double             x[3], y[3], a[3], b[3], c[3] ;
        double             angle = axis->order ? 2*M_PI/axis->order : 1.0 ;
        double             a_sin = sin( angle ) ;
        double             a_cos = cos( angle ) ;
        double             dot ;
        int                i ;

if( DIMENSION != 3 ){
    fprintf( stderr, "Catastrophe in rotate_atom!\n" ) ;
    exit( EXIT_FAILURE ) ;
    }
for( i = 0 ; i < 3 ; i++ )
    x[i] = from->x[i] - axis->distance * axis->normal[i] ;
for( i = 0, dot = 0 ; i < 3 ; i++ )
    dot += x[i] * axis->direction[i] ;
for( i = 0 ; i < 3 ; i++ )
    a[i] = axis->direction[i] * dot ;
for( i = 0 ; i < 3 ; i++ )
    b[i] = x[i] - a[i] ;
c[0] = b[1]*axis->direction[2] - b[2]*axis->direction[1] ;
c[1] = b[2]*axis->direction[0] - b[0]*axis->direction[2] ;
c[2] = b[0]*axis->direction[1] - b[1]*axis->direction[0] ;
for( i = 0 ; i < 3 ; i++ )
    y[i] = a[i] + b[i]*a_cos + c[i]*a_sin ;
for( i = 0 ; i < 3 ; i++ )
    to->x[i] = y[i] + axis->distance * axis->normal[i] ;
to->type = from->type ;
}

SYMMETRY_ELEMENT *
init_ultimate_axis(void)
{
        SYMMETRY_ELEMENT * axis = alloc_symmetry_element() ;
        double             dir[ DIMENSION ], rel[ DIMENSION ] ;
        double             s ;
        int                i, k ;

if( verbose > 0 ) printf( "Trying infinity axis\n" ) ;
StatTotal++ ;
axis->transform_atom = rotate_atom ;
axis->order          = 0 ;
axis->nparam         = 7 ;
for( k = 0 ; k < DIMENSION ; k++ )
    dir[k] = 0 ;
for( i = 0 ; i < AtomsCount ; i++ ){
    for( k = 0, s = 0 ; k < DIMENSION ; k++ ){
        rel[k] = Atoms[i].x[k] - CenterOfSomething[k] ;
        s     += rel[k]*dir[k] ;
        }
    if( s >= 0 )
         for( k = 0 ; k < DIMENSION ; k++ )
             dir[k] += rel[k] ;
    else for( k = 0 ; k < DIMENSION ; k++ )
             dir[k] -= rel[k] ;
    }
for( k = 0, s = 0 ; k < DIMENSION ; k++ )
    s += pow2( dir[k] ) ;
s = sqrt(s) ;
if( s > 0 )
     for( k = 0 ; k < DIMENSION ; k++ )
         dir[k] /= s ;
else dir[0] = 1 ;
for( k = 0 ; k < DIMENSION ; k++ )
    axis->direction[k] = dir[k] ;
for( k = 0, s = 0 ; k < DIMENSION ; k++ )
    s += pow2( CenterOfSomething[k] ) ;
s = sqrt(s) ;
if( s > 0 )
    for( k = 0 ; k < DIMENSION ; k++ )
        axis->normal[k] = CenterOfSomething[k]/s ;
else {
    for( k = 1 ; k < DIMENSION ; k++ )
        axis->normal[k] = 0 ;
    axis->normal[0] = 1 ;
    }
axis->distance = s ;
for( k = 0 ; k < AtomsCount ; k++ )
    axis->transform[k] = k ;
if( refine_symmetry_element( axis, 0 ) < 0 ){
    if( verbose > 0 ) printf( "    refinement failed for the infinity axis\n" ) ;
    destroy_symmetry_element( axis ) ;
    return NULL ;
    }
return axis ;
}


SYMMETRY_ELEMENT *
init_c2_axis( int i, int j, double support[ DIMENSION ] )
{
        SYMMETRY_ELEMENT * axis ;
        int                k ;
        double             ris, rjs ;
        double             r, center[ DIMENSION ] ;

if( verbose > 0 ) 
    printf( "Trying c2 axis for the pair (%d,%d) with the support (%g,%g,%g)\n", 
             i, j, support[0], support[1], support[2] ) ;
StatTotal++ ;
/* First, do a quick sanity check */
for( k = 0, ris = rjs = 0 ; k < DIMENSION ; k++ ){
    ris += pow2( Atoms[i].x[k] - support[k] ) ;
    rjs += pow2( Atoms[j].x[k] - support[k] ) ;
    }
ris = sqrt( ris ) ;
rjs = sqrt( rjs ) ;
if( fabs( ris - rjs ) > TolerancePrimary ){
    StatEarly++ ;
    if( verbose > 0 ) printf( "    Support can't actually define a rotation axis\n" ) ;
    return NULL ;
    }
axis                 = alloc_symmetry_element() ;
axis->transform_atom = rotate_atom ;
axis->order          = 2 ;
axis->nparam         = 7 ;
for( k = 0, r = 0 ; k < DIMENSION ; k++ )
    r += CenterOfSomething[k]*CenterOfSomething[k] ;
r = sqrt(r) ;
if( r > 0 ){
    for( k = 0 ; k < DIMENSION ; k++ )
        axis->normal[k] = CenterOfSomething[k]/r ;
    }
else {
    axis->normal[0] = 1 ;
    for( k = 1 ; k < DIMENSION ; k++ )
        axis->normal[k] = 0 ;
    }
axis->distance = r ;
for( k = 0, r = 0 ; k < DIMENSION ; k++ ){
    center[k] = ( Atoms[i].x[k] + Atoms[j].x[k] ) / 2 - support[k] ;
    r        += center[k]*center[k] ;
    }
r = sqrt(r) ;
if( r <= TolerancePrimary ){ /* c2 is underdefined, let's do something special */
    if( MolecularPlane != NULL ){
        if( verbose > 0 ) printf( "    c2 is underdefined, but there is a molecular plane\n" ) ;
        for( k = 0 ; k < DIMENSION ; k++ )
            axis->direction[k] = MolecularPlane->normal[k] ;
        }
    else {
        if( verbose > 0 ) printf( "    c2 is underdefined, trying random direction\n" ) ;
        for( k = 0 ; k < DIMENSION ; k++ )
            center[k] = Atoms[i].x[k] - Atoms[j].x[k] ;
        if( fabs( center[2] ) + fabs( center[1] ) > ToleranceSame ){
            axis->direction[0] =  0 ;
            axis->direction[1] =  center[2] ;
            axis->direction[2] = -center[1] ;
            }
        else {
            axis->direction[0] = -center[2] ;
            axis->direction[1] =  0 ;
            axis->direction[2] =  center[0] ;
            }
        for( k = 0, r = 0 ; k < DIMENSION ; k++ )
            r += axis->direction[k] * axis->direction[k] ;
        r = sqrt(r) ;
        for( k = 0 ; k < DIMENSION ; k++ )
            axis->direction[k] /= r ;
        }
    }
else { /* direction is Ok, renormalize it */
    for( k = 0 ; k < DIMENSION ; k++ )
        axis->direction[k] = center[k]/r ;
    }
if( refine_symmetry_element( axis, 1 ) < 0 ){
    if( verbose > 0 ) printf( "    refinement failed for the c2 axis\n" ) ;
    destroy_symmetry_element( axis ) ;
    return NULL ;
    }
return axis ;
}

SYMMETRY_ELEMENT *
init_axis_parameters( double a[3], double b[3], double c[3] )
{
        SYMMETRY_ELEMENT * axis ;
        int                i, order, sign ;
        double             ra, rb, rc, rab, rbc, rac, r ;
        double             angle ;

ra = rb = rc = rab = rbc = rac = 0 ;
for( i = 0 ; i < DIMENSION ; i++ ){
    ra  += a[i]*a[i] ;
    rb  += b[i]*b[i] ;
    rc  += c[i]*c[i] ;
    }
ra = sqrt(ra) ; rb  = sqrt(rb) ; rc  = sqrt(rc) ;
if( fabs( ra - rb ) > TolerancePrimary || fabs( ra - rc ) > TolerancePrimary || fabs( rb - rc ) > TolerancePrimary ){
    StatEarly++ ;
    if( verbose > 0 ) printf( "    points are not on a sphere\n" ) ;
    return NULL ;
    }
for( i = 0 ; i < DIMENSION ; i++ ){
    rab += (a[i]-b[i])*(a[i]-b[i]) ;
    rac += (a[i]-c[i])*(a[i]-c[i]) ;
    rbc += (c[i]-b[i])*(c[i]-b[i]) ;
    }
rab = sqrt(rab) ;
rac = sqrt(rac) ;
rbc = sqrt(rbc) ;
if( fabs( rab - rbc ) > TolerancePrimary ){
    StatEarly++ ;
    if( verbose > 0 ) printf( "    points can't be rotation-equivalent\n" ) ;
    return NULL ;
    }
if( rab <= ToleranceSame || rbc <= ToleranceSame || rac <= ToleranceSame ){
    StatEarly++ ;
    if( verbose > 0 ) printf( "    rotation is underdefined by these points\n" ) ;
    return NULL ;
    }
rab   = (rab+rbc)/2 ;
angle = M_PI - 2*asin( rac/(2*rab) ) ;
if( verbose > 1 ) printf( "    rotation angle is %f\n", angle ) ;
if( fabs(angle) <= M_PI/(MaxAxisOrder+1) ){
    StatEarly++ ;
    if( verbose > 0 ) printf( "    atoms are too close to a straight line\n" ) ;
    return NULL ;
    }
order = floor( (2*M_PI)/angle + 0.5 ) ;
if( order <= 2 || order > MaxAxisOrder ){
    StatEarly++ ;
    if( verbose > 0 ) printf( "    rotation axis order (%d) is not from 3 to %d\n", order, MaxAxisOrder ) ;
    return NULL ;
    }
axis = alloc_symmetry_element() ;
axis->order          = order ;
axis->nparam         = 7 ;
for( i = 0, r = 0 ; i < DIMENSION ; i++ )
    r += CenterOfSomething[i]*CenterOfSomething[i] ;
r = sqrt(r) ;
if( r > 0 ){
    for( i = 0 ; i < DIMENSION ; i++ )
        axis->normal[i] = CenterOfSomething[i]/r ;
    }
else {
    axis->normal[0] = 1 ;
    for( i = 1 ; i < DIMENSION ; i++ )
        axis->normal[i] = 0 ;
    }
axis->distance = r ;
axis->direction[0] = (b[1]-a[1])*(c[2]-b[2]) - (b[2]-a[2])*(c[1]-b[1]) ;
axis->direction[1] = (b[2]-a[2])*(c[0]-b[0]) - (b[0]-a[0])*(c[2]-b[2]) ;
axis->direction[2] = (b[0]-a[0])*(c[1]-b[1]) - (b[1]-a[1])*(c[0]-b[0]) ;
/*
 *  Arbitrarily select axis direction so that first non-zero component
 *  or the direction is positive.
 */
sign = 0 ;
if( axis->direction[0] <= 0 )
    if( axis->direction[0] < 0 )
         sign = 1 ;
    else if( axis->direction[1] <= 0 )
             if( axis->direction[1] < 0 )
                  sign = 1 ;
             else if( axis->direction[2] < 0 )
                      sign = 1 ;
if( sign )
    for( i = 0 ; i < DIMENSION ; i++ )
        axis->direction[i] = -axis->direction[i] ;
for( i = 0, r = 0 ; i < DIMENSION ; i++ )
    r += axis->direction[i]*axis->direction[i] ;
r = sqrt(r) ;
for( i = 0 ; i < DIMENSION ; i++ )
    axis->direction[i] /= r ;
if( verbose > 1 ){
    printf( "    axis origin is at (%g,%g,%g)\n", 
        axis->normal[0]*axis->distance, axis->normal[1]*axis->distance, axis->normal[2]*axis->distance ) ;
    printf( "    axis is in the direction (%g,%g,%g)\n", axis->direction[0], axis->direction[1], axis->direction[2] ) ;
    }
return axis ;
}

SYMMETRY_ELEMENT *
init_higher_axis( int ia, int ib, int ic )
{
        SYMMETRY_ELEMENT * axis ;
        double             a[ DIMENSION ], b[ DIMENSION ], c[ DIMENSION ] ;
        int                i ;

if( verbose > 0 ) printf( "Trying cn axis for the triplet (%d,%d,%d)\n", ia, ib, ic ) ;
StatTotal++ ;
/* Do a quick check of geometry validity */
for( i = 0 ; i < DIMENSION ; i++ ){
    a[i] = Atoms[ia].x[i] - CenterOfSomething[i] ;
    b[i] = Atoms[ib].x[i] - CenterOfSomething[i] ;
    c[i] = Atoms[ic].x[i] - CenterOfSomething[i] ;
    }
if( ( axis = init_axis_parameters( a, b, c ) ) == NULL ){
    if( verbose > 0 ) printf( "    no coherrent axis is defined by the points\n" ) ;
    return NULL ;
    }
axis->transform_atom = rotate_atom ;
if( refine_symmetry_element( axis, 1 ) < 0 ){
    if( verbose > 0 ) printf( "    refinement failed for the c%d axis\n", axis->order ) ;
    destroy_symmetry_element( axis ) ;
    return NULL ;
    }
return axis ;
}

/*
 *   Improper axes-specific routines.
 *   These are obtained by slight modifications of normal rotation
 *       routines.
 */
void
rotate_reflect_atom( SYMMETRY_ELEMENT *axis, ATOM1 *from, ATOM1 *to )
{
        double             x[3], y[3], a[3], b[3], c[3] ;
        double             angle = 2*M_PI/axis->order ;
        double             a_sin = sin( angle ) ;
        double             a_cos = cos( angle ) ;
        double             dot ;
        int                i ;

if( DIMENSION != 3 ){
    fprintf( stderr, "Catastrophe in rotate_reflect_atom!\n" ) ;
    exit( EXIT_FAILURE ) ;
    }
for( i = 0 ; i < 3 ; i++ )
    x[i] = from->x[i] - axis->distance * axis->normal[i] ;
for( i = 0, dot = 0 ; i < 3 ; i++ )
    dot += x[i] * axis->direction[i] ;
for( i = 0 ; i < 3 ; i++ )
    a[i] = axis->direction[i] * dot ;
for( i = 0 ; i < 3 ; i++ )
    b[i] = x[i] - a[i] ;
c[0] = b[1]*axis->direction[2] - b[2]*axis->direction[1] ;
c[1] = b[2]*axis->direction[0] - b[0]*axis->direction[2] ;
c[2] = b[0]*axis->direction[1] - b[1]*axis->direction[0] ;
for( i = 0 ; i < 3 ; i++ )
    y[i] = -a[i] + b[i]*a_cos + c[i]*a_sin ;
for( i = 0 ; i < 3 ; i++ )
    to->x[i] = y[i] + axis->distance * axis->normal[i] ;
to->type = from->type ;
}

SYMMETRY_ELEMENT *
init_improper_axis( int ia, int ib, int ic )
{
        SYMMETRY_ELEMENT * axis ;
        double             a[ DIMENSION ], b[ DIMENSION ], c[ DIMENSION ] ;
        double             centerpoint[ DIMENSION ] ;
        double             r ;
        int                i ;

if( verbose > 0 ) printf( "Trying sn axis for the triplet (%d,%d,%d)\n", ia, ib, ic ) ;
StatTotal++ ;
/* First, reduce the problem to Cn case */
for( i = 0 ; i < DIMENSION ; i++ ){
    a[i] = Atoms[ia].x[i] - CenterOfSomething[i] ;
    b[i] = Atoms[ib].x[i] - CenterOfSomething[i] ;
    c[i] = Atoms[ic].x[i] - CenterOfSomething[i] ;
    }
for( i = 0, r = 0 ; i < DIMENSION ; i++ ){
    centerpoint[i] = a[i] + c[i] + 2*b[i] ;
    r             += centerpoint[i]*centerpoint[i] ;
    }
r = sqrt(r) ;
if( r <= ToleranceSame ){
    StatEarly++ ;
    if( verbose > 0 ) printf( "    atoms can not define improper axis of the order more than 2\n" ) ;
    return NULL ;
    }
for( i = 0 ; i < DIMENSION ; i++ )
    centerpoint[i] /= r ;
for( i = 0, r = 0 ; i < DIMENSION ; i++ )
    r += centerpoint[i] * b[i] ;
for( i = 0 ; i < DIMENSION ; i++ )
    b[i] = 2*r*centerpoint[i] - b[i] ;
/* Do a quick check of geometry validity */
if( ( axis = init_axis_parameters( a, b, c ) ) == NULL ){
    if( verbose > 0 ) printf( "    no coherrent improper axis is defined by the points\n" ) ;
    return NULL ;
    }
axis->transform_atom = rotate_reflect_atom ;
if( refine_symmetry_element( axis, 1 ) < 0 ){
    if( verbose > 0 ) printf( "    refinement failed for the s%d axis\n", axis->order ) ;
    destroy_symmetry_element( axis ) ;
    return NULL ;
    }
return axis ;
}

/*
 *   Control routines
 */

void
find_center_of_something( void )
{
        int                i, j ;
        double             coord_sum[ DIMENSION ] ;
        double             r ;

for( j = 0 ; j < DIMENSION ; j++ )
    coord_sum[j] = 0 ;
for( i = 0 ; i < AtomsCount ; i++ ){
    for( j = 0 ; j < DIMENSION ; j++ )
        coord_sum[j] += Atoms[i].x[j] ;
    }
for( j = 0 ; j < DIMENSION ; j++ )
    CenterOfSomething[j] = coord_sum[j]/AtomsCount ;
if( verbose > 0 )
    printf( "Center of something is at %15.10f, %15.10f, %15.10f\n", 
            CenterOfSomething[0], CenterOfSomething[1], CenterOfSomething[2] ) ;
DistanceFromCenter = (double *) calloc( AtomsCount, sizeof( double ) ) ;
if( DistanceFromCenter == NULL ){
    fprintf( stderr, "Unable to allocate array for the distances\n" ) ;
    exit( EXIT_FAILURE ) ;
    }
for( i = 0 ; i < AtomsCount ; i++ ){
    for( j = 0, r = 0 ; j < DIMENSION ; j++ )
        r += pow2( Atoms[i].x[j] - CenterOfSomething[j] ) ;
    DistanceFromCenter[i] = r ;
    }
}

void
find_planes(void)
{
        int                i, j ;
        SYMMETRY_ELEMENT * plane ;

plane = init_ultimate_plane() ;
if( plane != NULL ){
    MolecularPlane = plane ;
    PlanesCount++ ;
    Planes = (SYMMETRY_ELEMENT **) realloc( Planes, sizeof( SYMMETRY_ELEMENT* ) * PlanesCount ) ;
    if( Planes == NULL ){
        perror( "Out of memory in find_planes" ) ;
        exit( EXIT_FAILURE ) ;
        }
    Planes[ PlanesCount - 1 ] = plane ;
    }
for( i = 1 ; i < AtomsCount ; i++ ){
    for( j = 0 ; j < i ; j++ ){
        if( Atoms[i].type != Atoms[j].type )
            continue ;
        if( ( plane = init_mirror_plane( i, j ) ) != NULL ){
            PlanesCount++ ;
            Planes = (SYMMETRY_ELEMENT **) realloc( Planes, sizeof( SYMMETRY_ELEMENT* ) * PlanesCount ) ;
            if( Planes == NULL ){
                perror( "Out of memory in find_planes" ) ;
                exit( EXIT_FAILURE ) ;
                }
            Planes[ PlanesCount - 1 ] = plane ;
            }
        }
    }
}

void
find_inversion_centers(void)
{
        SYMMETRY_ELEMENT * center ;

if( ( center = init_inversion_center() ) != NULL ){
    InversionCenters = (SYMMETRY_ELEMENT **) calloc( 1, sizeof( SYMMETRY_ELEMENT* ) ) ;
    InversionCenters[0]   = center ;
    InversionCentersCount = 1 ;
    }
}

void
find_infinity_axis(void)
{
        SYMMETRY_ELEMENT * axis ;

if( ( axis = init_ultimate_axis() ) != NULL ){
    NormalAxesCount++ ;
    NormalAxes = (SYMMETRY_ELEMENT **) realloc( NormalAxes, sizeof( SYMMETRY_ELEMENT* ) * NormalAxesCount ) ;
    if( NormalAxes == NULL ){
        perror( "Out of memory in find_infinity_axes()" ) ;
        exit( EXIT_FAILURE ) ;
        }
    NormalAxes[ NormalAxesCount - 1 ] = axis ;
    }
}

void
find_c2_axes(void)
{
        int                i, j, k, l, m ;
        double             center[ DIMENSION ] ;
        double *           distances = calloc( AtomsCount, sizeof( double ) ) ;
        double             r ;
        SYMMETRY_ELEMENT * axis ;

if( distances == NULL ){
    fprintf( stderr, "Out of memory in find_c2_axes()\n" ) ;
    exit( EXIT_FAILURE ) ;
    }
for( i = 1 ; i < AtomsCount ; i++ ){
    for( j = 0 ; j < i ; j++ ){
        if( Atoms[i].type != Atoms[j].type )
            continue ;
        if( fabs( DistanceFromCenter[i] - DistanceFromCenter[j] ) > TolerancePrimary )
            continue ; /* A very cheap, but quite effective check */
        /*
         *   First, let's try to get it cheap and use CenterOfSomething
         */
        for( k = 0, r = 0 ; k < DIMENSION ; k++ ){
            center[k] = ( Atoms[i].x[k] + Atoms[j].x[k] ) / 2 ;
            r        += pow2( center[k] - CenterOfSomething[k] ) ;
            }
        r = sqrt(r) ;
        if( r > 5*TolerancePrimary ){ /* It's Ok to use CenterOfSomething */
            if( ( axis = init_c2_axis( i, j, CenterOfSomething ) ) != NULL ){
                NormalAxesCount++ ;
                NormalAxes = (SYMMETRY_ELEMENT **) realloc( NormalAxes, sizeof( SYMMETRY_ELEMENT* ) * NormalAxesCount ) ;
                if( NormalAxes == NULL ){
                    perror( "Out of memory in find_c2_axes" ) ;
                    exit( EXIT_FAILURE ) ;
                    }
                NormalAxes[ NormalAxesCount - 1 ] = axis ;
                }
            continue ;
            }
        /*
         *  Now, C2 axis can either pass through an atom, or through the
         *  middle of the other pair.
         */
        for( k = 0 ; k < AtomsCount ; k++ ){
            if( ( axis = init_c2_axis( i, j, Atoms[k].x ) ) != NULL ){
                NormalAxesCount++ ;
                NormalAxes = (SYMMETRY_ELEMENT **) realloc( NormalAxes, sizeof( SYMMETRY_ELEMENT* ) * NormalAxesCount ) ;
                if( NormalAxes == NULL ){
                    perror( "Out of memory in find_c2_axes" ) ;
                    exit( EXIT_FAILURE ) ;
                    }
                NormalAxes[ NormalAxesCount - 1 ] = axis ;
                }
            }
        /*
         *  Prepare data for an additional pre-screening check
         */
        for( k = 0 ; k < AtomsCount ; k++ ){
            for( l = 0, r = 0 ; l < DIMENSION ; l++ )
                r += pow2( Atoms[k].x[l] - center[l] ) ;
            distances[k] = sqrt(r) ;
            }
        for( k = 0 ; k < AtomsCount ; k++ ){
            for( l = 0 ; l < AtomsCount ; l++ ){
                if( Atoms[k].type != Atoms[l].type )
                    continue ;
                if( fabs( DistanceFromCenter[k] - DistanceFromCenter[l] ) > TolerancePrimary ||
                    fabs( distances[k] - distances[l] ) > TolerancePrimary )
                        continue ; /* We really need this one to run reasonably fast! */
                for( m = 0 ; m < DIMENSION ; m++ )
                    center[m] = ( Atoms[k].x[m] + Atoms[l].x[m] ) / 2 ;
                if( ( axis = init_c2_axis( i, j, center ) ) != NULL ){
                    NormalAxesCount++ ;
                    NormalAxes = (SYMMETRY_ELEMENT **) realloc( NormalAxes, sizeof( SYMMETRY_ELEMENT* ) * NormalAxesCount ) ;
                    if( NormalAxes == NULL ){
                        perror( "Out of memory in find_c2_axes" ) ;
                        exit( EXIT_FAILURE ) ;
                        }
                    NormalAxes[ NormalAxesCount - 1 ] = axis ;
                    }
                }
            }
        }
    }
free( distances ) ;
}

void
find_higher_axes(void)
{
        int                i, j, k ;
        SYMMETRY_ELEMENT * axis ;

for( i = 0 ; i < AtomsCount ; i++ ){
    for( j = i + 1 ; j < AtomsCount ; j++ ){
        if( Atoms[i].type != Atoms[j].type )
            continue ;
        if( fabs( DistanceFromCenter[i] - DistanceFromCenter[j] ) > TolerancePrimary )
            continue ; /* A very cheap, but quite effective check */
        for( k = 0 ; k < AtomsCount ; k++ ){
            if( Atoms[i].type != Atoms[k].type )
                continue ;
            if( ( fabs( DistanceFromCenter[i] - DistanceFromCenter[k] ) > TolerancePrimary ) ||
                ( fabs( DistanceFromCenter[j] - DistanceFromCenter[k] ) > TolerancePrimary ) )
                    continue ;
            if( ( axis = init_higher_axis( i, j, k ) ) != NULL ){
                NormalAxesCount++ ;
                NormalAxes = (SYMMETRY_ELEMENT **) realloc( NormalAxes, sizeof( SYMMETRY_ELEMENT* ) * NormalAxesCount ) ;
                if( NormalAxes == NULL ){
                    perror( "Out of memory in find_higher_axes" ) ;
                    exit( EXIT_FAILURE ) ;
                    }
                NormalAxes[ NormalAxesCount - 1 ] = axis ;
                }
            }
        }
    }
}

void
find_improper_axes(void)
{
        int                i, j, k ;
        SYMMETRY_ELEMENT * axis ;

for( i = 0 ; i < AtomsCount ; i++ ){
    for( j = i + 1 ; j < AtomsCount ; j++ ){
        for( k = 0 ; k < AtomsCount ; k++ ){
            if( ( axis = init_improper_axis( i, j, k ) ) != NULL ){
                ImproperAxesCount++ ;
                ImproperAxes = (SYMMETRY_ELEMENT **) realloc( ImproperAxes, sizeof( SYMMETRY_ELEMENT* ) * ImproperAxesCount ) ;
                if( ImproperAxes == NULL ){
                    perror( "Out of memory in find_higher_axes" ) ;
                    exit( EXIT_FAILURE ) ;
                    }
                ImproperAxes[ ImproperAxesCount - 1 ] = axis ;
                }
            }
        }
    }
}

void
report_planes( void )
{
        int           i ;

if( PlanesCount == 0 )
    printf( "There are no planes of symmetry in the molecule\n" ) ;
else {
    if( PlanesCount == 1 )
         printf( "There is a plane of symmetry in the molecule\n" ) ;
    else printf( "There are %d planes of symmetry in the molecule\n", PlanesCount ) ;
    printf( "    Rank Residual          Direction of the normal           Distance\n" ) ;
    for( i = 0 ; i < PlanesCount ; i++ ){
        printf( "%3d%4d %8.4e ", i+1, Planes[i]->rank, Planes[i]->maxdev ) ;
        printf( "(%11.8f,%11.8f,%11.8f) ", Planes[i]->normal[0], Planes[i]->normal[1], Planes[i]->normal[2] ) ;
        printf( "%14.8f\n", Planes[i]->distance ) ;
        }
    }
}

void
report_inversion_centers( void )
{
if( InversionCentersCount == 0 )
     printf( "There is no inversion center in the molecule\n" ) ;
else {
    printf( "There is an inversion center in the molecule\n" ) ;
    printf( " Rank Residual                      Position\n" ) ;
    printf( "%4d %8.4e ", InversionCenters[0]->rank, InversionCenters[0]->maxdev ) ;
    printf( "(%14.8f,%14.8f,%14.8f)\n",
        InversionCenters[0]->distance * InversionCenters[0]->normal[0],
        InversionCenters[0]->distance * InversionCenters[0]->normal[1],
        InversionCenters[0]->distance * InversionCenters[0]->normal[2] ) ;
    }
}

void
report_axes( void )
{
        int           i ;

if( NormalAxesCount == 0 )
    printf( "There are no normal axes in the molecule\n" ) ;
else {
    if( NormalAxesCount == 1 )
         printf( "There is a normal axis in the molecule\n" ) ;
    else printf( "There are %d normal axes in the molecule\n", NormalAxesCount ) ;
    printf( "    Rank Residual  Order         Direction of the axis                         Supporting point\n" ) ;
    for( i = 0 ; i < NormalAxesCount ; i++ ){
        printf( "%3d%4d %8.4e ", i+1, NormalAxes[i]->rank, NormalAxes[i]->maxdev ) ;
        if( NormalAxes[i]->order == 0 )
             printf( "Inf " ) ;
        else printf( "%3d ", NormalAxes[i]->order ) ;
        printf( "(%11.8f,%11.8f,%11.8f) ", 
            NormalAxes[i]->direction[0], NormalAxes[i]->direction[1], NormalAxes[i]->direction[2] ) ;
        printf( "(%14.8f,%14.8f,%14.8f)\n", 
            NormalAxes[0]->distance * NormalAxes[0]->normal[0],
            NormalAxes[0]->distance * NormalAxes[0]->normal[1],
            NormalAxes[0]->distance * NormalAxes[0]->normal[2] ) ;
        }
    }
}

void
report_improper_axes( void )
{
        int           i ;

if( ImproperAxesCount == 0 )
    printf( "There are no improper axes in the molecule\n" ) ;
else {
    if( ImproperAxesCount == 1 )
         printf( "There is an improper axis in the molecule\n" ) ;
    else printf( "There are %d improper axes in the molecule\n", ImproperAxesCount ) ;
    printf( "    Rank Residual  Order         Direction of the axis                         Supporting point\n" ) ;
    for( i = 0 ; i < ImproperAxesCount ; i++ ){
        printf( "%3d%4d %8.4e ", i+1, ImproperAxes[i]->rank, ImproperAxes[i]->maxdev ) ;
        if( ImproperAxes[i]->order == 0 )
             printf( "Inf " ) ;
        else printf( "%3d ", ImproperAxes[i]->order ) ;
        printf( "(%11.8f,%11.8f,%11.8f) ", 
            ImproperAxes[i]->direction[0], ImproperAxes[i]->direction[1], ImproperAxes[i]->direction[2] ) ;
        printf( "(%14.8f,%14.8f,%14.8f)\n", 
            ImproperAxes[0]->distance * ImproperAxes[0]->normal[0],
            ImproperAxes[0]->distance * ImproperAxes[0]->normal[1],
            ImproperAxes[0]->distance * ImproperAxes[0]->normal[2] ) ;
        }
    }
}

/*
 *  General symmetry handling
 */
void
report_and_reset_counters( void )
{
printf( "  %10ld candidates examined\n"
        "  %10ld removed early\n"
        "  %10ld removed during initial mating stage\n"
        "  %10ld removed as duplicates\n"
        "  %10ld removed because of the wrong transformation order\n"
        "  %10ld removed after unsuccessful optimization\n"
        "  %10ld accepted\n",
    StatTotal, StatEarly, StatPairs, StatDups, StatOrder, StatOpt, StatAccept ) ;
StatTotal = StatEarly = StatPairs = StatDups = StatOrder = StatOpt = StatAccept = 0 ;
}

void
find_symmetry_elements( void )
{
find_center_of_something() ;
if( verbose > -1 ){
    printf( "Looking for the inversion center\n" ) ;
    }
find_inversion_centers() ;
if( verbose > -1 ){
    report_and_reset_counters() ;
    printf( "Looking for the planes of symmetry\n" ) ;
    }
find_planes() ;
if( verbose > -1 ){
    report_and_reset_counters() ;
    printf( "Looking for infinity axis\n" ) ;
    }
find_infinity_axis() ;
if( verbose > -1 ){
    report_and_reset_counters() ;
    printf( "Looking for C2 axes\n" ) ;
    }
find_c2_axes() ;
if( verbose > -1 ){
    report_and_reset_counters() ;
    printf( "Looking for higher axes\n" ) ;
    }
find_higher_axes() ;
if( verbose > -1 ){
    report_and_reset_counters() ;
    printf( "Looking for the improper axes\n" ) ;
    }
find_improper_axes() ;
if( verbose > -1 ){
    report_and_reset_counters() ;
    }
}

int
compare_axes( const void *a, const void *b )
{
        SYMMETRY_ELEMENT * axis_a = *(SYMMETRY_ELEMENT**) a ;
        SYMMETRY_ELEMENT * axis_b = *(SYMMETRY_ELEMENT**) b ;
        int                i, order_a, order_b ;

order_a = axis_a->order ; if( order_a == 0 ) order_a = 10000 ;
order_b = axis_b->order ; if( order_b == 0 ) order_b = 10000 ;
if( ( i = order_b - order_a ) != 0 ) return i ;
if( axis_a->maxdev > axis_b->maxdev ) return -1 ;
if( axis_a->maxdev < axis_b->maxdev ) return  1 ;
return 0 ;
}

void
sort_symmetry_elements( void )
{
if( PlanesCount > 1 ){
    qsort( Planes, PlanesCount, sizeof( SYMMETRY_ELEMENT * ), compare_axes ) ;
    }
if( NormalAxesCount > 1 ){
    qsort( NormalAxes, NormalAxesCount, sizeof( SYMMETRY_ELEMENT * ), compare_axes ) ;
    }
if( ImproperAxesCount > 1 ){
    qsort( ImproperAxes, ImproperAxesCount, sizeof( SYMMETRY_ELEMENT * ), compare_axes ) ;
    }
}

void
report_symmetry_elements_verbose( void )
{
report_inversion_centers() ;
report_axes() ;
report_improper_axes() ;
report_planes() ;
}

void
summarize_symmetry_elements( void )
{
        int          i ;

NormalAxesCounts   = (int*) calloc( MaxAxisOrder+1, sizeof( int ) ) ;
ImproperAxesCounts = (int*) calloc( MaxAxisOrder+1, sizeof( int ) ) ;
for( i = 0 ; i < NormalAxesCount ; i++ )
    NormalAxesCounts[ NormalAxes[i]->order ]++ ;
for( i = 0 ; i < ImproperAxesCount ; i++ )
    ImproperAxesCounts[ ImproperAxes[i]->order ]++ ;
}

void
report_symmetry_elements_brief( void )
{
        int          i ;
        char *       symmetry_code = calloc( 1, 10*(PlanesCount+NormalAxesCount+ImproperAxesCount+InversionCentersCount+2) ) ;
        char         buf[ 100 ] ;

if( symmetry_code == NULL ){
    if (verbose > -3)
        fprintf( stderr, "Unable to allocate memory for symmetry ID code in report_symmetry_elements_brief()\n" ) ;
    exit( EXIT_FAILURE ) ;
    }
if (PlanesCount + NormalAxesCount + ImproperAxesCount + InversionCentersCount == 0)
    {
        if (verbose > -3)
            printf("Molecule has no symmetry elements\n");
    }
else {
    if (verbose > -3)
        printf( "Molecule has the following symmetry elements: " ) ;
    if( InversionCentersCount > 0 ) strcat( symmetry_code, "(i) " ) ;
    if( NormalAxesCounts[0] == 1 )
         strcat( symmetry_code, "(Cinf) " ) ;
    if( NormalAxesCounts[0] >  1 ) {
        sprintf( buf, "%d*(Cinf) ", NormalAxesCounts[0] ) ;
        strcat( symmetry_code, buf ) ;
        }
    for( i = MaxAxisOrder ; i >= 2 ; i-- ){
        if( NormalAxesCounts[i] == 1 ){ sprintf( buf, "(C%d) ", i ) ; strcat( symmetry_code, buf ) ; }
        if( NormalAxesCounts[i] >  1 ){ sprintf( buf, "%d*(C%d) ", NormalAxesCounts[i], i ) ; strcat( symmetry_code, buf ) ; }
        }
    for( i = MaxAxisOrder ; i >= 2 ; i-- ){
        if( ImproperAxesCounts[i] == 1 ){ sprintf( buf, "(S%d) ", i ) ; strcat( symmetry_code, buf ) ; }
        if( ImproperAxesCounts[i] >  1 ){ sprintf( buf, "%d*(S%d) ", ImproperAxesCounts[i], i ) ; strcat( symmetry_code, buf ) ; }
        }
    if( PlanesCount == 1 ) strcat( symmetry_code, "(sigma) " ) ;
    if( PlanesCount >  1 ){ sprintf( buf, "%d*(sigma) ", PlanesCount ) ; strcat( symmetry_code, buf ) ; }
    if (verbose > -3)
        printf( "%s\n", symmetry_code ) ;
    }
SymmetryCode = symmetry_code ;
}

void
identify_point_group( void )
{
        int            i ;
        int            last_matching = -1 ;
        int            matching_count = 0 ;

for( i = 0 ; i < PointGroupsCount ; i++ ){
    if( strcmp( SymmetryCode, PointGroups[i].symmetry_code ) == 0 ){
        if( PointGroups[i].check() == 1 ){
            last_matching = i ;
            matching_count++ ;
            }
        else {
            if( verbose > -2 ){
                printf( "It looks very much like %s, but it is not since %s\n", 
                    PointGroups[i].group_name, PointGroupRejectionReason ) ;
                }
            }
        }
    }
if( matching_count == 0 ){
    printf( "These symmetry elements match no point group I know of. Sorry.\n" ) ;
    }
if( matching_count >  1 ){
    printf( "These symmetry elements match more than one group I know of.\n"
            "SOMETHING IS VERY WRONG\n" ) ;
    printf( "Matching groups are:\n" ) ;
    for( i = 0 ; i < PointGroupsCount ; i++ ){
        if( ( strcmp( SymmetryCode, PointGroups[i].symmetry_code ) == 0 ) && ( PointGroups[i].check() == 1 ) ){
            printf( "    %s\n", PointGroups[i].group_name ) ;
            }
        }
    }
if( matching_count == 1 ){
    printf( "It seems to be the %s point group\n", PointGroups[last_matching].group_name ) ;
    }
}

void
identify_point_group1( char *PG )
{
        int            i ;
        int            last_matching = -1 ;
        int            matching_count = 0 ;

for( i = 0 ; i < PointGroupsCount ; i++ ){
    if( strcmp( SymmetryCode, PointGroups[i].symmetry_code ) == 0 ){
        if( PointGroups[i].check() == 1 ){
            last_matching = i ;
            matching_count++ ;
            }
        else {
            if( verbose > -2 ){
                printf( "It looks very much like %s, but it is not since %s\n", 
                    PointGroups[i].group_name, PointGroupRejectionReason ) ;
                }
            }
        }
    }
if( matching_count == 0 ){
    if (verbose > -2)
        printf( "These symmetry elements match no point group I know of. Sorry.\n" ) ;
    }
if( matching_count >  1 ){
    if (verbose > -2)
        printf( "These symmetry elements match more than one group I know of.\n"
            "SOMETHING IS VERY WRONG\n" ) ;
    if (verbose > -2)
        printf( "Matching groups are:\n" ) ;
    for( i = 0 ; i < PointGroupsCount ; i++ ){
        if( ( strcmp( SymmetryCode, PointGroups[i].symmetry_code ) == 0 ) && ( PointGroups[i].check() == 1 ) ){
            if (verbose > -2)
                printf( "    %s\n", PointGroups[i].group_name ) ;
            }
        }
    }
if( matching_count == 1 ){
    if (verbose > -2)
        printf( "It seems to be the %s point group\n", PointGroups[last_matching].group_name ) ;
    }
strcpy(PG, PointGroups[last_matching].group_name);

}



/*
 *  Input/Output
 */

int
read_coordinates( FILE *in )
{
        int                 i ;

if( fscanf( in, "%d", &AtomsCount ) != 1 ){
    fprintf( stderr, "Error reading atom count\n" ) ;
    return -1 ;
    }
if( verbose > 0 ) printf( "Atoms count = %d\n", AtomsCount ) ;
Atoms = calloc( AtomsCount, sizeof( ATOM1 ) ) ;
if( Atoms == NULL ){
    fprintf( stderr, "Out of memory for atoms coordinates\n" ) ;
    return -1 ;
    }
for( i = 0 ; i < AtomsCount ; i++ ){
    if( fscanf( in, "%d %lg %lg %lg\n", &Atoms[i].type, &Atoms[i].x[0], &Atoms[i].x[1], &Atoms[i].x[2] ) != 4 ){
        fprintf( stderr, "Error reading description of the atom %d\n", i ) ;
        return -1 ;
        }
    }
return 0 ;
}

/*************************************************************************
 * Begin SGE code
 *************************************************************************/

/*************************************************************************
 * Normalize vector and align it to the axis
 * if its direction differs not more than the tolerance
 *************************************************************************/
static void normalize (double *norm)
{
  int		k;
  double	dist;

  for (k=0; k<DIMENSION; k++)
  {
    if (fabs(norm[k])   <= ToleranceFinal) norm[k] = 0;
    if (fabs(norm[k]-1) <= ToleranceFinal) norm[k] = 1;
    if (fabs(norm[k]+1) <= ToleranceFinal) norm[k] = -1;
  }
  dist = 0;
  for (k=0; k<DIMENSION; k++) dist += norm[k]*norm[k];
  dist = sqrt (dist);
  if (dist > ToleranceFinal) for (k=0; k<DIMENSION; k++) norm[k] /= dist;
  else for (k=0; k<DIMENSION; k++) norm[k] = 0;
}

/*************************************************************************
 * Translate coordinates of all atoms and symmetry elements
 * to the distance of dist in the direction of the norm vector
 *************************************************************************/
static void translate_coords (double dist, double *norm)
{
  int i, k;
  
  /* For all atoms... */
  for (i=0; i<AtomsCount; i++)
    for (k=0; k<DIMENSION; k++) Atoms[i].x[k] -= dist*norm[k];

  /* For all planes... */
  for (i=0; i<PlanesCount; i++)
  {
    for (k=0; k<DIMENSION; k++)			/* vector subtraction */
      Planes[i]->distance -= dist*Planes[i]->normal[k]*norm[k];
    if (fabs(Planes[i]->distance) <= ToleranceFinal) Planes[i]->distance = 0;
    if (Planes[i]->distance >= 0) continue;
    Planes[i]->distance = -Planes[i]->distance;	/* make distance positive */
    for (k=0; k<DIMENSION; k++) Planes[i]->normal[k] = -Planes[i]->normal[k];
  }

  /* For all normal axes... */
  for (i=0; i<NormalAxesCount; i++)
  {
    for (k=0; k<DIMENSION; k++)		/* transform normal into coords */
    {
      NormalAxes[i]->normal[k] =
	NormalAxes[i]->normal[k]*NormalAxes[i]->distance-dist*norm[k];
      if (fabs(NormalAxes[i]->normal[k]) <= ToleranceFinal)
	NormalAxes[i]->normal[k] = 0;
    }
    NormalAxes[i]->distance = 0;
    for (k=0; k<DIMENSION; k++)
      NormalAxes[i]->distance +=
	NormalAxes[i]->normal[k]*NormalAxes[i]->normal[k];
    NormalAxes[i]->distance = sqrt(NormalAxes[i]->distance);
    if (NormalAxes[i]->distance > ToleranceFinal)	/* normalize it back */
    {
      for (k=0; k<DIMENSION; k++)
	NormalAxes[i]->normal[k] /= NormalAxes[i]->distance;
      normalize (NormalAxes[i]->normal);
    }
    else				/* dont normalize zero distance */
    {
      NormalAxes[i]->distance = 0;
      for (k=0; k<DIMENSION; k++) NormalAxes[i]->normal[k] = 0;
    }
  }

  /* For all improper axes... */
  for (i=0; i<ImproperAxesCount; i++)
  {
    for (k=0; k<DIMENSION; k++)		/* transform normal into coords */
    {
      ImproperAxes[i]->normal[k] =
	ImproperAxes[i]->normal[k]*ImproperAxes[i]->distance-dist*norm[k];
      if (fabs(ImproperAxes[i]->normal[k]) <= ToleranceFinal)
	ImproperAxes[i]->normal[k] = 0;
    }
    ImproperAxes[i]->distance = 0;
    for (k=0; k<DIMENSION; k++)
      ImproperAxes[i]->distance +=
	ImproperAxes[i]->normal[k]*ImproperAxes[i]->normal[k];
    ImproperAxes[i]->distance = sqrt(ImproperAxes[i]->distance);
    if (ImproperAxes[i]->distance > ToleranceFinal)	/* normalize it back */
    {
      for (k=0; k<DIMENSION; k++)
	ImproperAxes[i]->normal[k] /= ImproperAxes[i]->distance;
      normalize (ImproperAxes[i]->normal);
    }
    else				/* dont normalize zero distance */
    {
      ImproperAxes[i]->distance = 0;
      for (k=0; k<DIMENSION; k++) ImproperAxes[i]->normal[k] = 0;
    }
  }

  /* For all inversion centers... */
  for (i=0; i<InversionCentersCount; i++)
  {
    for (k=0; k<DIMENSION; k++)		/* transform normal into coords */
    {
      InversionCenters[i]->normal[k] =
	InversionCenters[i]->normal[k]*InversionCenters[i]->distance-
	dist*norm[k];
      if (fabs(InversionCenters[i]->normal[k]) <= ToleranceFinal)
	InversionCenters[i]->normal[k] = 0;
    }
    InversionCenters[i]->distance = 0;
    for (k=0; k<DIMENSION; k++)
      InversionCenters[i]->distance +=
	InversionCenters[i]->normal[k]*InversionCenters[i]->normal[k];
    InversionCenters[i]->distance = sqrt(InversionCenters[i]->distance);
    if (InversionCenters[i]->distance > ToleranceFinal)	/* normalize it back */
    {
      for (k=0; k<DIMENSION; k++)
	InversionCenters[i]->normal[k] /= InversionCenters[i]->distance;
      normalize (InversionCenters[i]->normal);
    }
    else				/* dont normalize zero distance */
    {
      InversionCenters[i]->distance = 0;
      for (k=0; k<DIMENSION; k++) InversionCenters[i]->normal[k] = 0;
    }
  }
}

/*************************************************************************
 * Rotate coordinates of all atoms and symmetry elements so that
 * the direction of the norm vector be aligned to the specified axis
 *************************************************************************/
static void rotate_coords (int axis, double *norm)
{
  int		i, k;
  double	dist, x, y, z, norm1[DIMENSION], norm2[DIMENSION];

  /*
   * Make orthogonal coordinate system for future rotation
   * norm1[] will be perpendicular to norm[] and to the specified axis
   * norm2[] will be perpendicular to norm[] and to norm1[]
   */
  switch (axis)
  {
  case 0:								/* X */
    norm1[0] = 0;
    norm1[1] = norm[2];
    norm1[2] = -norm[1];
    norm2[0] = -norm[1]*norm[1]-norm[2]*norm[2];
    norm2[1] = norm[0]*norm[1];
    norm2[2] = norm[0]*norm[2];
    break;
  case 1:								/* Y */
    norm1[0] = -norm[2];
    norm1[1] = 0;
    norm1[2] = norm[0];
    norm2[0] = norm[1]*norm[0];
    norm2[1] = -norm[2]*norm[2]-norm[0]*norm[0];
    norm2[2] = norm[1]*norm[2];
    break;
  case 2:								/* Z */
    norm1[0] = norm[1];
    norm1[1] = -norm[0];
    norm1[2] = 0;
    norm2[0] = norm[2]*norm[0];
    norm2[1] = norm[2]*norm[1];
    norm2[2] = -norm[0]*norm[0]-norm[1]*norm[1];
    break;
  default:							/* ??? */
    fprintf (stderr, "Catastrophe in rotate_coords !\n");
    exit (EXIT_FAILURE);
  }

  /* Normalize norm1 and norm2... */
  dist = 0;
  for (k=0; k<DIMENSION; k++) dist += norm1[k]*norm1[k];
  dist = sqrt (dist);
  for (k=0; k<DIMENSION; k++) norm1[k] /= dist;
  dist = 0;
  for (k=0; k<DIMENSION; k++) dist += norm2[k]*norm2[k];
  dist = sqrt (dist);
  for (k=0; k<DIMENSION; k++) norm2[k] /= dist;

  /*
   * Make target coordinate system for future rotation swapping X, Y, Z
   * so that (norm,norm1,norm2) be as near to (X,Y,Z) as possible
   */
  switch (axis)
  {
  case 0:						/* already OK */
    break;
  case 1:						/* change to (Z,X,Y) */
    for (k=0; k<DIMENSION; k++)
    {
      x = norm[k];
      y = norm1[k];
      z = norm2[k];
      norm[k]  = z;
      norm1[k] = x;
      norm2[k] = y;
    }
    break;
  case 2:						/* change to (Y,Z,X) */
    for (k=0; k<DIMENSION; k++)
    {
      x = norm[k];
      y = norm1[k];
      z = norm2[k];
      norm[k]  = y;
      norm1[k] = z;
      norm2[k] = x;
    }
    break;
  default:							/* ??? */
    fprintf (stderr, "Catastrophe in rotate_coords !\n");
    exit (EXIT_FAILURE);
  }

  /* Now rotate all atoms... */
  for (i=0; i<AtomsCount; i++)
  {
    x = y = z = 0;
    for (k=0; k<DIMENSION; k++)
    {
      x += Atoms[i].x[k]*norm[k];
      y += Atoms[i].x[k]*norm1[k];
      z += Atoms[i].x[k]*norm2[k];
    }
    Atoms[i].x[0] = x;
    Atoms[i].x[1] = y;
    Atoms[i].x[2] = z;
  }

  /* Rotate all planes... */
  for (i=0; i<PlanesCount; i++)
  {
    x = y = z = 0;
    for (k=0; k<DIMENSION; k++)
    {
      x += Planes[i]->normal[k]*norm[k];
      y += Planes[i]->normal[k]*norm1[k];
      z += Planes[i]->normal[k]*norm2[k];
    }
    Planes[i]->normal[0] = x;
    Planes[i]->normal[1] = y;
    Planes[i]->normal[2] = z;
    normalize (Planes[i]->normal);
  }

  /* All normal axes... */
  for (i=0; i<NormalAxesCount; i++)
  {
    x = y = z = 0;
    for (k=0; k<DIMENSION; k++)
    {
      x += NormalAxes[i]->normal[k]*norm[k];
      y += NormalAxes[i]->normal[k]*norm1[k];
      z += NormalAxes[i]->normal[k]*norm2[k];
    }
    NormalAxes[i]->normal[0] = x;
    NormalAxes[i]->normal[1] = y;
    NormalAxes[i]->normal[2] = z;
    normalize (NormalAxes[i]->normal);
    x = y = z = 0;
    for (k=0; k<DIMENSION; k++)
    {
      x += NormalAxes[i]->direction[k]*norm[k];
      y += NormalAxes[i]->direction[k]*norm1[k];
      z += NormalAxes[i]->direction[k]*norm2[k];
    }
    NormalAxes[i]->direction[0] = x;
    NormalAxes[i]->direction[1] = y;
    NormalAxes[i]->direction[2] = z;
    normalize (NormalAxes[i]->direction);
  }

  /* All improper axes... */
  for (i=0; i<ImproperAxesCount; i++)
  {
    x = y = z = 0;
    for (k=0; k<DIMENSION; k++)
    {
      x += ImproperAxes[i]->normal[k]*norm[k];
      y += ImproperAxes[i]->normal[k]*norm1[k];
      z += ImproperAxes[i]->normal[k]*norm2[k];
    }
    ImproperAxes[i]->normal[0] = x;
    ImproperAxes[i]->normal[1] = y;
    ImproperAxes[i]->normal[2] = z;
    normalize (ImproperAxes[i]->normal);
    x = y = z = 0;
    for (k=0; k<DIMENSION; k++)
    {
      x += ImproperAxes[i]->direction[k]*norm[k];
      y += ImproperAxes[i]->direction[k]*norm1[k];
      z += ImproperAxes[i]->direction[k]*norm2[k];
    }
    ImproperAxes[i]->direction[0] = x;
    ImproperAxes[i]->direction[1] = y;
    ImproperAxes[i]->direction[2] = z;
    normalize (ImproperAxes[i]->direction);
  }

  /* All inversion centers... */
  for (i=0; i<InversionCentersCount; i++)
  {
    x = y = z = 0;
    for (k=0; k<DIMENSION; k++)
    {
      x += InversionCenters[i]->normal[k]*norm[k];
      y += InversionCenters[i]->normal[k]*norm1[k];
      z += InversionCenters[i]->normal[k]*norm2[k];
    }
    InversionCenters[i]->normal[0] = x;
    InversionCenters[i]->normal[1] = y;
    InversionCenters[i]->normal[2] = z;
    normalize (InversionCenters[i]->normal);
  }
}

/*************************************************************************
 * Swap coordinates of all atoms and symmetry elements
 * so that the specified coordinate become the coordinate Z
 *************************************************************************/
static void swap_coords (int coord)
{
  int		i;
  double	tmp;
  
  if (coord == 2) return;		/* wanted coordinate is Z already */

  /* For all atoms... */
  for (i=0; i<AtomsCount; i++)
  {
    tmp = Atoms[i].x[2];
    Atoms[i].x[2] = Atoms[i].x[coord];
    Atoms[i].x[coord] = tmp;
  }

  /* For all planes... */
  for (i=0; i<PlanesCount; i++)
  {
    tmp = Planes[i]->normal[2];
    Planes[i]->normal[2] = Planes[i]->normal[coord];
    Planes[i]->normal[coord] = tmp;
  }

  /* For all normal axes... */
  for (i=0; i<NormalAxesCount; i++)
  {
    tmp = NormalAxes[i]->normal[2];
    NormalAxes[i]->normal[2] = NormalAxes[i]->normal[coord];
    NormalAxes[i]->normal[coord] = tmp;
    tmp = NormalAxes[i]->direction[2];
    NormalAxes[i]->direction[2] = NormalAxes[i]->direction[coord];
    NormalAxes[i]->direction[coord] = tmp;
  }

  /* For all improper axes... */
  for (i=0; i<ImproperAxesCount; i++)
  {
    tmp = ImproperAxes[i]->normal[2];
    ImproperAxes[i]->normal[2] = ImproperAxes[i]->normal[coord];
    ImproperAxes[i]->normal[coord] = tmp;
    tmp = ImproperAxes[i]->direction[2];
    ImproperAxes[i]->direction[2] = ImproperAxes[i]->direction[coord];
    ImproperAxes[i]->direction[coord] = tmp;
  }

  /* For all inversion centers... */
  for (i=0; i<InversionCentersCount; i++)
  {
    tmp = InversionCenters[i]->normal[2];
    InversionCenters[i]->normal[2] = InversionCenters[i]->normal[coord];
    InversionCenters[i]->normal[coord] = tmp;
  }
}

/*************************************************************************
 * Swap coordinates of all atoms and symmetry elements
 * so that the specified coordinate become the coordinate Y
 *************************************************************************/
static void swap_y_coords (int coord)
{
  int		i;
  double	tmp;
  
  if (coord == 1) return;		/* wanted coordinate is Y already */

  /* For all atoms... */
  for (i=0; i<AtomsCount; i++)
  {
    tmp = Atoms[i].x[1];
    Atoms[i].x[1] = Atoms[i].x[coord];
    Atoms[i].x[coord] = tmp;
  }

  /* For all planes... */
  for (i=0; i<PlanesCount; i++)
  {
    tmp = Planes[i]->normal[1];
    Planes[i]->normal[1] = Planes[i]->normal[coord];
    Planes[i]->normal[coord] = tmp;
  }

  /* For all normal axes... */
  for (i=0; i<NormalAxesCount; i++)
  {
    tmp = NormalAxes[i]->normal[1];
    NormalAxes[i]->normal[1] = NormalAxes[i]->normal[coord];
    NormalAxes[i]->normal[coord] = tmp;
    tmp = NormalAxes[i]->direction[1];
    NormalAxes[i]->direction[1] = NormalAxes[i]->direction[coord];
    NormalAxes[i]->direction[coord] = tmp;
  }

  /* For all improper axes... */
  for (i=0; i<ImproperAxesCount; i++)
  {
    tmp = ImproperAxes[i]->normal[1];
    ImproperAxes[i]->normal[1] = ImproperAxes[i]->normal[coord];
    ImproperAxes[i]->normal[coord] = tmp;
    tmp = ImproperAxes[i]->direction[1];
    ImproperAxes[i]->direction[1] = ImproperAxes[i]->direction[coord];
    ImproperAxes[i]->direction[coord] = tmp;
  }

  /* For all inversion centers... */
  for (i=0; i<InversionCentersCount; i++)
  {
    tmp = InversionCenters[i]->normal[1];
    InversionCenters[i]->normal[1] = InversionCenters[i]->normal[coord];
    InversionCenters[i]->normal[coord] = tmp;
  }
}

/*************************************************************************
 * Swap coordinates of all atoms and symmetry elements
 * so that the specified coordinate become the coordinate X
 *************************************************************************/
static void swap_x_coords (int coord)
{
  int		i;
  double	tmp;
  
  if (coord == 0) return;		/* wanted coordinate is X already */

  /* For all atoms... */
  for (i=0; i<AtomsCount; i++)
  {
    tmp = Atoms[i].x[0];
    Atoms[i].x[0] = Atoms[i].x[coord];
    Atoms[i].x[coord] = tmp;
  }

  /* For all planes... */
  for (i=0; i<PlanesCount; i++)
  {
    tmp = Planes[i]->normal[0];
    Planes[i]->normal[0] = Planes[i]->normal[coord];
    Planes[i]->normal[coord] = tmp;
  }

  /* For all normal axes... */
  for (i=0; i<NormalAxesCount; i++)
  {
    tmp = NormalAxes[i]->normal[0];
    NormalAxes[i]->normal[0] = NormalAxes[i]->normal[coord];
    NormalAxes[i]->normal[coord] = tmp;
    tmp = NormalAxes[i]->direction[0];
    NormalAxes[i]->direction[0] = NormalAxes[i]->direction[coord];
    NormalAxes[i]->direction[coord] = tmp;
  }

  /* For all improper axes... */
  for (i=0; i<ImproperAxesCount; i++)
  {
    tmp = ImproperAxes[i]->normal[0];
    ImproperAxes[i]->normal[0] = ImproperAxes[i]->normal[coord];
    ImproperAxes[i]->normal[coord] = tmp;
    tmp = ImproperAxes[i]->direction[0];
    ImproperAxes[i]->direction[0] = ImproperAxes[i]->direction[coord];
    ImproperAxes[i]->direction[coord] = tmp;
  }

  /* For all inversion centers... */
  for (i=0; i<InversionCentersCount; i++)
  {
    tmp = InversionCenters[i]->normal[0];
    InversionCenters[i]->normal[0] = InversionCenters[i]->normal[coord];
    InversionCenters[i]->normal[coord] = tmp;
  }
}

/*************************************************************************
 * Reorient symmetry planes so that they go through the coordinate origin
 * and be perpendicular to any of the coordinate axes if possible
 * Before translation/rotation it is checked if no previously aligned
 * symmetry plane can become misaligned in the result of the operation
 *************************************************************************/
static void reorient_planes (void)
{
  int		i, j, k, i_one, i_max, n_one, n_zero;
  double	dist, dmax, norm[DIMENSION], norm1[DIMENSION];
  
  for (i=0; i<PlanesCount; i++)				/* for all planes */
  {
    /* Firstly translate the plane... */
    do						/* dummy loop to avoid goto */
    {
      if (fabs(Planes[i]->distance) <= ToleranceFinal) Planes[i]->distance = 0;
      dist = Planes[i]->distance;
      normalize (Planes[i]->normal);
      for (k=0; k<DIMENSION; k++) norm[k] = Planes[i]->normal[k];
      if (dist <= ToleranceFinal) break;	/* already acceptable */

      dmax = 0;
      for (j=0; j<i; j++)		/* check previously aligned planes */
      {
	normalize (Planes[j]->normal);
	for (k=0; k<DIMENSION; k++) norm1[k] = Planes[j]->normal[k];
	n_one = n_zero = 0;
	for (k=0; k<DIMENSION; k++)
	{
	  if (fabs(fabs(norm1[k])-1) <= ToleranceFinal) n_one  ++;
	  if (fabs(norm1[k])         <= ToleranceFinal) n_zero ++;
	}
	if (n_one!=1 || n_zero!=2) continue;	/* this plane unaligned */
	for (k=0; k<DIMENSION; k++) dmax += norm[k]*norm1[k];
	if (fabs(dmax) > ToleranceFinal &&
	    OrdUniq*Planes[j]->rank <= OrdUniq*Planes[i]->rank)
	  break;				/* aligned non-perpendicular */
	dmax = 0;		/* allowed to translate perpendicular plane */
      }
      if (fabs(dmax) > ToleranceFinal) break;	/* there is something bad */

      translate_coords (dist, norm);		/* translation allowed */
    }
    while (0);				/* perform dummy loop only once */

    /* Now rotate the plane... */
    normalize (Planes[i]->normal);
    for (k=0; k<DIMENSION; k++) norm[k] = Planes[i]->normal[k];
    i_one = i_max  = -1;
    n_one = n_zero = 0;
    dmax  = -1;
    for (k=0; k<DIMENSION; k++)
    {
      if (fabs(norm[k]) > dmax)		/* look for the nearest coordinate */
      {
	i_max = k;
	dmax  = fabs(norm[k]);
      }
      if (fabs(fabs(norm[k])-1) <= ToleranceFinal)
      {
	i_one = k;
	n_one ++;
      }
      if (fabs(norm[k]) <= ToleranceFinal) n_zero ++;
    }
    if (n_one==1 && n_zero==2) continue;	/* already acceptable */

    dmax = 0;
    for (j=0; j<i; j++)			/* check previously aligned planes */
    {
      normalize (Planes[j]->normal);
      for (k=0; k<DIMENSION; k++) norm1[k] = Planes[j]->normal[k];
      n_one = n_zero = 0;
      for (k=0; k<DIMENSION; k++)
      {
	if (fabs(fabs(norm1[k])-1) <= ToleranceFinal) n_one  ++;
	if (fabs(norm1[k])         <= ToleranceFinal) n_zero ++;
      }
      if (n_one!=1 || n_zero!=2) continue;	/* this plane unaligned */
      for (k=0; k<DIMENSION; k++) dmax += norm[k]*norm1[k];
      if (fabs(dmax) > ToleranceFinal &&
	  OrdUniq*Planes[j]->rank <= OrdUniq*Planes[i]->rank)
	break;					/* aligned non-perpendicular */
      dmax = 0;			/* allowed to rotate perpendicular plane */
    }
    if (fabs(dmax) > ToleranceFinal) continue;	/* there is something bad */

    rotate_coords (i_max, norm);			/* rotation allowed */
  }							/* for all planes */
}

/*************************************************************************
 * Reorient normal axes so that they go through the coordinate origin
 * and be parallel to any of the coordinate axes if possible
 * Similar to reorient_planes()
 * Before translation/rotation it is checked if no previously aligned
 * symmetry plane or normal axis can become misaligned
 * in the result of the operation
 * As an exception abelian axis (order 2 or inf) may misalign non-abelian axes
 *************************************************************************/
static void reorient_axes (int abel)
{
  int		i, j, k, i_one, i_max, n_one, n_zero;
  double	dist, dmax, norm[DIMENSION], norm1[DIMENSION];
  
  for (i=0; i<NormalAxesCount; i++)			/* for all axes */
  {
    /* Firstly translate the axis... */
    do						/* dummy loop to avoid goto */
    {
      if (fabs(NormalAxes[i]->distance) <= ToleranceFinal)
      {
	NormalAxes[i]->distance = 0;
	for (k=0; k<DIMENSION; k++) NormalAxes[i]->normal[k] = 0;
      }
      dist = NormalAxes[i]->distance;
      normalize (NormalAxes[i]->normal);
      for (k=0; k<DIMENSION; k++) norm[k] = NormalAxes[i]->normal[k];
      if (dist <= ToleranceFinal) break;	/* already acceptable */

      dmax = 0;
      for (j=0; j<PlanesCount; j++)	/* check previously aligned planes */
      {
	normalize (Planes[j]->normal);
	for (k=0; k<DIMENSION; k++) norm1[k] = Planes[j]->normal[k];
	n_one = n_zero = 0;
	for (k=0; k<DIMENSION; k++)
	{
	  if (fabs(fabs(norm1[k])-1) <= ToleranceFinal) n_one  ++;
	  if (fabs(norm1[k])         <= ToleranceFinal) n_zero ++;
	}
	if (n_one!=1 || n_zero!=2) continue;	/* this plane unaligned */
	for (k=0; k<DIMENSION; k++) dmax += norm[k]*norm1[k];
	if (fabs(dmax) > ToleranceFinal &&
	    (abel ||
	     (NormalAxes[i]->order && NormalAxes[i]->order <= OrdAxis*2)))
	  break;				/* aligned non-perpendicular */
	dmax = 0;		/* allowed to translate perpendicular plane */
      }
      if (fabs(dmax) > ToleranceFinal) break;	/* there is something bad */

      dmax = 1;
      for (j=0; j<i; j++)		/* check previously aligned axes */
      {
	normalize (NormalAxes[j]->direction);
	for (k=0; k<DIMENSION; k++) norm1[k] = NormalAxes[j]->direction[k];
	n_one = n_zero = 0;
	for (k=0; k<DIMENSION; k++)
	{
	  if (fabs(fabs(norm1[k])-1) <= ToleranceFinal) n_one  ++;
	  if (fabs(norm1[k])         <= ToleranceFinal) n_zero ++;
	}
	if (n_one!=1 || n_zero!=2) continue;	/* this axis unaligned */
	dmax = 0;
	for (k=0; k<DIMENSION; k++) dmax += norm[k]*norm1[k];
	if (fabs(fabs(dmax)-1) > ToleranceFinal &&
	    ((abel &&
	      ((NormalAxes[i]->order && NormalAxes[i]->order != 2) ||
	       !NormalAxes[j]->order || NormalAxes[j]->order == 2) &&
	      OrdUniq*NormalAxes[j]->rank <= OrdUniq*NormalAxes[i]->rank) ||
	     (!abel &&
	      (!NormalAxes[j]->order ||
	       OrdAxis*NormalAxes[j]->order > OrdAxis*NormalAxes[i]->order ||
	       (NormalAxes[j]->order == NormalAxes[i]->order &&
		OrdUniq*NormalAxes[j]->rank <= OrdUniq*NormalAxes[i]->rank)))))
	  break;			/* aligned non-parallel abelian */
	dmax = 1;	/* allowed to translate parallel or non-abelian axis */
      }
      if (fabs(fabs(dmax)-1) > ToleranceFinal) break;/*there is something bad*/

      translate_coords (dist, norm);		/* translation allowed */
    }
    while (0);				/* perform dummy loop only once */

    /* Now rotate the axis... */
    normalize (NormalAxes[i]->direction);
    for (k=0; k<DIMENSION; k++) norm[k] = NormalAxes[i]->direction[k];
    i_one = i_max  = -1;
    n_one = n_zero = 0;
    dmax  = -1;
    for (k=0; k<DIMENSION; k++)
    {
      if (fabs(norm[k]) > dmax)		/* look for the nearest coordinate */
      {
	i_max = k;
	dmax  = fabs(norm[k]);
      }
      if (fabs(fabs(norm[k])-1) <= ToleranceFinal)
      {
	i_one = k;
	n_one ++;
      }
      if (fabs(norm[k]) <= ToleranceFinal) n_zero ++;
    }
    if (n_one==1 && n_zero==2) continue;	/* already acceptable */

    dmax = 0;
    for (j=0; j<PlanesCount; j++)	/* check previously aligned planes */
    {
      normalize (Planes[j]->normal);
      for (k=0; k<DIMENSION; k++) norm1[k] = Planes[j]->normal[k];
      n_one = n_zero = 0;
      for (k=0; k<DIMENSION; k++)
      {
	if (fabs(fabs(norm1[k])-1) <= ToleranceFinal) n_one  ++;
	if (fabs(norm1[k])         <= ToleranceFinal) n_zero ++;
      }
      if (n_one!=1 || n_zero!=2) continue;	/* this plane unaligned */
      for (k=0; k<DIMENSION; k++) dmax += norm[k]*norm1[k];
      if (fabs(dmax) > ToleranceFinal &&
	  (abel ||
	   (NormalAxes[i]->order && NormalAxes[i]->order <= OrdAxis*2)))
	break;					/* aligned non-perpendicular */
      dmax = 0;			/* allowed to rotate perpendicular plane */
    }
    if (fabs(dmax) > ToleranceFinal) continue;	/* there is something bad */

    dmax = 0;
    for (j=0; j<i; j++)			/* check previously aligned axes */
    {
      normalize (NormalAxes[j]->direction);
      for (k=0; k<DIMENSION; k++) norm1[k] = NormalAxes[j]->direction[k];
      n_one = n_zero = 0;
      for (k=0; k<DIMENSION; k++)
      {
	if (fabs(fabs(norm1[k])-1) <= ToleranceFinal) n_one  ++;
	if (fabs(norm1[k])         <= ToleranceFinal) n_zero ++;
      }
      if (n_one!=1 || n_zero!=2) continue;	/* this axis unaligned */
      for (k=0; k<DIMENSION; k++) dmax += norm[k]*norm1[k];
      if (fabs(dmax) > ToleranceFinal &&
	  ((abel &&
	    ((NormalAxes[i]->order && NormalAxes[i]->order != 2) ||
	     !NormalAxes[j]->order || NormalAxes[j]->order == 2) &&
	    OrdUniq*NormalAxes[j]->rank <= OrdUniq*NormalAxes[i]->rank) ||
	   (!abel &&
	    (!NormalAxes[j]->order ||
	     OrdAxis*NormalAxes[j]->order > OrdAxis*NormalAxes[i]->order ||
	     (NormalAxes[j]->order == NormalAxes[i]->order &&
	      OrdUniq*NormalAxes[j]->rank <= OrdUniq*NormalAxes[i]->rank)))))
	break;				/* aligned non-perpendicular abelian */
      dmax = 0;	/* allowed to rotate perpendicular or non-abelian axis */
    }
    if (fabs(dmax) > ToleranceFinal) continue;	/* there is something bad */

    rotate_coords (i_max, norm);			/* rotation allowed */
  }							/* for all axes */
}

/*************************************************************************
 * Reorient improper axes so that they go through the coordinate origin
 * and be parallel to any of the coordinate axes if possible
 * Similar to reorient_axes()
 * Before translation/rotation it is checked if no previously aligned
 * symmetry plane or axis can become misaligned in the result of the operation
 *************************************************************************/
static void reorient_improper_axes (void)
{
  int		i, j, k, i_one, i_max, n_one, n_zero;
  double	dist, dmax, norm[DIMENSION], norm1[DIMENSION];
  
  for (i=0; i<ImproperAxesCount; i++)			/* for all axes */
  {
    /* Firstly translate the axis... */
    do						/* dummy loop to avoid goto */
    {
      if (fabs(ImproperAxes[i]->distance) <= ToleranceFinal)
      {
	ImproperAxes[i]->distance = 0;
	for (k=0; k<DIMENSION; k++) ImproperAxes[i]->normal[k] = 0;
      }
      dist = ImproperAxes[i]->distance;
      normalize (ImproperAxes[i]->normal);
      for (k=0; k<DIMENSION; k++) norm[k] = ImproperAxes[i]->normal[k];
      if (dist <= ToleranceFinal) break;	/* already acceptable */

      dmax = 0;
      for (j=0; j<PlanesCount; j++)	/* check previously aligned planes */
      {
	normalize (Planes[j]->normal);
	for (k=0; k<DIMENSION; k++) norm1[k] = Planes[j]->normal[k];
	n_one = n_zero = 0;
	for (k=0; k<DIMENSION; k++)
	{
	  if (fabs(fabs(norm1[k])-1) <= ToleranceFinal) n_one  ++;
	  if (fabs(norm1[k])         <= ToleranceFinal) n_zero ++;
	}
	if (n_one!=1 || n_zero!=2) continue;	/* this plane unaligned */
	for (k=0; k<DIMENSION; k++) dmax += norm[k]*norm1[k];
	if (fabs(dmax) > ToleranceFinal) break;	/* aligned non-perpendicular */
	dmax = 0;		/* allowed to translate perpendicular plane */
      }
      if (fabs(dmax) > ToleranceFinal) break;	/* there is something bad */

      dmax = 1;
      for (j=0; j<NormalAxesCount; j++)	/* check previously aligned axes */
      {
	normalize (NormalAxes[j]->direction);
	for (k=0; k<DIMENSION; k++) norm1[k] = NormalAxes[j]->direction[k];
	n_one = n_zero = 0;
	for (k=0; k<DIMENSION; k++)
	{
	  if (fabs(fabs(norm1[k])-1) <= ToleranceFinal) n_one  ++;
	  if (fabs(norm1[k])         <= ToleranceFinal) n_zero ++;
	}
	if (n_one!=1 || n_zero!=2) continue;	/* this axis unaligned */
	dmax = 0;
	for (k=0; k<DIMENSION; k++) dmax += norm[k]*norm1[k];
	if (fabs(fabs(dmax)-1) > ToleranceFinal) break;/*aligned non-parallel*/
	dmax = 1;		/* allowed to translate parallel axis */
      }
      if (fabs(fabs(dmax)-1) > ToleranceFinal) break;	/* something bad */

      dmax = 1;
      for (j=0; j<i; j++)	/* check previously aligned improper axes */
      {
	normalize (ImproperAxes[j]->direction);
	for (k=0; k<DIMENSION; k++) norm1[k] = ImproperAxes[j]->direction[k];
	n_one = n_zero = 0;
	for (k=0; k<DIMENSION; k++)
	{
	  if (fabs(fabs(norm1[k])-1) <= ToleranceFinal) n_one  ++;
	  if (fabs(norm1[k])         <= ToleranceFinal) n_zero ++;
	}
	if (n_one!=1 || n_zero!=2) continue; /* this improper axis unaligned */
	dmax = 0;
	for (k=0; k<DIMENSION; k++) dmax += norm[k]*norm1[k];
	if (fabs(fabs(dmax)-1) > ToleranceFinal &&
	    OrdUniq*ImproperAxes[j]->rank <= OrdUniq*ImproperAxes[i]->rank)
	  break;				/* aligned non-parallel */
	dmax = 1;	/* allowed to translate parallel improper axis */
      }
      if (fabs(fabs(dmax)-1) > ToleranceFinal) break;	/* something bad */

      translate_coords (dist, norm);		/* translation allowed */
    }
    while (0);				/* perform dummy loop only once */

    /* Now rotate the axis... */
    normalize (ImproperAxes[i]->direction);
    for (k=0; k<DIMENSION; k++) norm[k] = ImproperAxes[i]->direction[k];
    i_one = i_max  = -1;
    n_one = n_zero = 0;
    dmax  = -1;
    for (k=0; k<DIMENSION; k++)
    {
      if (fabs(norm[k]) > dmax)		/* look for the nearest coordinate */
      {
	i_max = k;
	dmax  = fabs(norm[k]);
      }
      if (fabs(fabs(norm[k])-1) <= ToleranceFinal)
      {
	i_one = k;
	n_one ++;
      }
      if (fabs(norm[k]) <= ToleranceFinal) n_zero ++;
    }
    if (n_one==1 && n_zero==2) continue;	/* already acceptable */

    dmax = 0;
    for (j=0; j<PlanesCount; j++)	/* check previously aligned planes */
    {
      normalize (Planes[j]->normal);
      for (k=0; k<DIMENSION; k++) norm1[k] = Planes[j]->normal[k];
      n_one = n_zero = 0;
      for (k=0; k<DIMENSION; k++)
      {
	if (fabs(fabs(norm1[k])-1) <= ToleranceFinal) n_one  ++;
	if (fabs(norm1[k])         <= ToleranceFinal) n_zero ++;
      }
      if (n_one!=1 || n_zero!=2) continue;	/* this plane unaligned */
      for (k=0; k<DIMENSION; k++) dmax += norm[k]*norm1[k];
      if (fabs(dmax) > ToleranceFinal) break;	/* aligned non-perpendicular */
      dmax = 0;			/* allowed to rotate perpendicular plane */
    }
    if (fabs(dmax) > ToleranceFinal) continue;	/* there is something bad */

    dmax = 0;
    for (j=0; j<NormalAxesCount; j++)	/* check previously aligned axes */
    {
      normalize (NormalAxes[j]->direction);
      for (k=0; k<DIMENSION; k++) norm1[k] = NormalAxes[j]->direction[k];
      n_one = n_zero = 0;
      for (k=0; k<DIMENSION; k++)
      {
	if (fabs(fabs(norm1[k])-1) <= ToleranceFinal) n_one  ++;
	if (fabs(norm1[k])         <= ToleranceFinal) n_zero ++;
      }
      if (n_one!=1 || n_zero!=2) continue;	/* this axis unaligned */
      for (k=0; k<DIMENSION; k++) dmax += norm[k]*norm1[k];
      if (fabs(dmax) > ToleranceFinal) break;	/* aligned non-perpendicular */
      dmax = 0;			/* allowed to rotate perpendicular axis */
    }
    if (fabs(dmax) > ToleranceFinal) continue;	/* there is something bad */

    dmax = 0;
    for (j=0; j<i; j++)		/* check previously aligned improper axes */
    {
      normalize (ImproperAxes[j]->direction);
      for (k=0; k<DIMENSION; k++) norm1[k] = ImproperAxes[j]->direction[k];
      n_one = n_zero = 0;
      for (k=0; k<DIMENSION; k++)
      {
	if (fabs(fabs(norm1[k])-1) <= ToleranceFinal) n_one  ++;
	if (fabs(norm1[k])         <= ToleranceFinal) n_zero ++;
      }
      if (n_one!=1 || n_zero!=2) continue;   /* this improper axis unaligned */
      for (k=0; k<DIMENSION; k++) dmax += norm[k]*norm1[k];
      if (fabs(dmax) > ToleranceFinal &&
	  OrdUniq*ImproperAxes[j]->rank <= OrdUniq*ImproperAxes[i]->rank)
	break;					/* aligned non-perpendicular */
      dmax = 0;		/* allowed to rotate perpendicular improper axis */
    }
    if (fabs(dmax) > ToleranceFinal) continue;	/* there is something bad */

    rotate_coords (i_max, norm);			/* rotation allowed */
  }							/* for all axes */
}

/*************************************************************************
 * Reorient inversion centers so that they be equal to the coordinate origin
 * Similar to reorient_axes() but without rotation
 * Before translation it is checked if no previously aligned symmetry plane
 * or normal axis can become misaligned in the result of the operation
 * As an exception it is allowed to misalign non-abelian axes
 *************************************************************************/
static void reorient_inversion_centers (void)
{
  int		i, j, k, n_one, n_zero;
  double	dist, dmax, norm[DIMENSION], norm1[DIMENSION];
  
  for (i=0; i<InversionCentersCount; i++)		/* for all centers */
  {
    if (fabs(InversionCenters[i]->distance) <= ToleranceFinal)
    {
      InversionCenters[i]->distance = 0;
      for (k=0; k<DIMENSION; k++) InversionCenters[i]->normal[k] = 0;
    }
    dist = InversionCenters[i]->distance;
    normalize (InversionCenters[i]->normal);
    for (k=0; k<DIMENSION; k++) norm[k] = InversionCenters[i]->normal[k];
    if (dist <= ToleranceFinal) continue;	/* already acceptable */

    dmax = 0;
    for (j=0; j<PlanesCount; j++)	/* check previously aligned planes */
    {
      normalize (Planes[j]->normal);
      for (k=0; k<DIMENSION; k++) norm1[k] = Planes[j]->normal[k];
      n_one = n_zero = 0;
      for (k=0; k<DIMENSION; k++)
      {
	if (fabs(fabs(norm1[k])-1) <= ToleranceFinal) n_one  ++;
	if (fabs(norm1[k])         <= ToleranceFinal) n_zero ++;
      }
      if (n_one!=1 || n_zero!=2) continue;	/* this plane unaligned */
      for (k=0; k<DIMENSION; k++) dmax += norm[k]*norm1[k];
      if (fabs(dmax) > ToleranceFinal) break;	/* aligned non-perpendicular */
      dmax = 0;			/* allowed to translate perpendicular plane */
    }
    if (fabs(dmax) > ToleranceFinal) continue;	/* there is something bad */

    dmax = 1;
    for (j=0; j<NormalAxesCount; j++)	/* check previously aligned axes */
    {
      normalize (NormalAxes[j]->direction);
      for (k=0; k<DIMENSION; k++) norm1[k] = NormalAxes[j]->direction[k];
      n_one = n_zero = 0;
      for (k=0; k<DIMENSION; k++)
      {
	if (fabs(fabs(norm1[k])-1) <= ToleranceFinal) n_one  ++;
	if (fabs(norm1[k])         <= ToleranceFinal) n_zero ++;
      }
      if (n_one!=1 || n_zero!=2) continue;	/* this axis unaligned */
      dmax = 0;
      for (k=0; k<DIMENSION; k++) dmax += norm[k]*norm1[k];
      if (fabs(fabs(dmax)-1) > ToleranceFinal &&
	  (!NormalAxes[j]->order || NormalAxes[j]->order == 2))
	break;				/* aligned non-parallel abelian */
      dmax = 1;		/* allowed to translate parallel or non-abelian axis */
    }
    if (fabs(fabs(dmax)-1) > ToleranceFinal &&
	(!NormalAxes[j]->order || NormalAxes[j]->order == 2))
      continue;					/* there is something bad */

    translate_coords (dist, norm);		/* translation allowed */
  }							/* for all centers */
}

/*************************************************************************
 * Select the proper coordinate axes for the master frame
 *************************************************************************/
static void canonize_axes (int abel)
{
  int		i, k, i_one, i_zero, n_one, n_zero, nplanes, naxes, maxis;
  int		cplanes[DIMENSION], caxes[DIMENSION];
  double	norm[DIMENSION];

  /* Align infinity axis (if any) along Z */
  naxes = 0;
  for (k=0; k<DIMENSION; k++) caxes[k] = 0;
  for (i=0; i<NormalAxesCount; i++)
  {
    if (NormalAxes[i]->order) continue;		/* not an infinity axis */
    normalize (NormalAxes[i]->direction);
    for (k=0; k<DIMENSION; k++) norm[k] = NormalAxes[i]->direction[k];
    i_one = -1;
    n_one = n_zero = 0;
    for (k=0; k<DIMENSION; k++)
    {
      if (fabs(fabs(norm[k])-1) <= ToleranceFinal)
      {
	i_one = k;
	n_one ++;
      }
      if (fabs(norm[k]) <= ToleranceFinal) n_zero ++;
    }
    if (n_one!=1 || n_zero!=2) continue;		/* not abelian axis */
    naxes ++;
    caxes[i_one] = 1;						/* count it */
    break;
  }
  if (naxes == 1)		/* align single infinity axis along Z */
  {
    i_one = 2;
    for (k=0; k<DIMENSION; k++) if (caxes[k]) i_one = k;       /* find coord */
    swap_coords (i_one);
    return;			/* Cinfv or Dinfh or Kh aka C2v or D2h */
  }

  /* Align principal axis (if any) along Z */
  if (! abel)
  {
    maxis = 0;
    for (i=0; i<NormalAxesCount; i++)
    {
      if (! NormalAxes[i]->order) continue;		/* an infinity axis */
      normalize (NormalAxes[i]->direction);
      for (k=0; k<DIMENSION; k++) norm[k] = NormalAxes[i]->direction[k];
      i_one = -1;
      n_one = n_zero = 0;
      for (k=0; k<DIMENSION; k++)
      {
	if (fabs(fabs(norm[k])-1) <= ToleranceFinal)
	{
	  i_one = k;
	  n_one ++;
	}
	if (fabs(norm[k]) <= ToleranceFinal) n_zero ++;
      }
      if (n_one!=1 || n_zero!=2) continue;		/* not aligned axis */
      if (maxis < NormalAxes[i]->order) maxis = NormalAxes[i]->order;
    }
    if (maxis > 2)					/* there exists one */
    {
      naxes = 0;
      for (k=0; k<DIMENSION; k++) caxes[k] = 0;
      for (i=0; i<NormalAxesCount; i++)
      {
	if (! NormalAxes[i]->order) continue;		/* an infinity axis */
	normalize (NormalAxes[i]->direction);
	for (k=0; k<DIMENSION; k++) norm[k] = NormalAxes[i]->direction[k];
	i_one = -1;
	n_one = n_zero = 0;
	for (k=0; k<DIMENSION; k++)
	{
	  if (fabs(fabs(norm[k])-1) <= ToleranceFinal)
	  {
	    i_one = k;
	    n_one ++;
	  }
	  if (fabs(norm[k]) <= ToleranceFinal) n_zero ++;
	}
	if (n_one!=1 || n_zero!=2)         continue;	/* not aligned axis */
	if (NormalAxes[i]->order != maxis) continue;   /* not principal axis */
	if (caxes[i_one])                  continue;	/* already counted */
	naxes ++;
	caxes[i_one] = 1;					/* count it */
      }
      if (naxes == 2)		/* align two order 2 axes perpendicular to Z */
      {
	i_zero = 2;
	for (k=0; k<DIMENSION; k++) if (! caxes[k]) i_zero = k;/* find coord */
	swap_coords (i_zero);
      }
      if (naxes == 1)		/* align single principal axis along Z */
      {
	i_one = 2;
	for (k=0; k<DIMENSION; k++) if (caxes[k]) i_one = k;   /* find coord */
	swap_coords (i_one);
      }
      /* Align perpendicular order 2 normal axes (if any) */
      naxes = 0;
      for (k=0; k<DIMENSION; k++) caxes[k] = 0;
      for (i=0; i<NormalAxesCount; i++)
      {
	if (NormalAxes[i]->order != 2) continue;      /* not an order 2 axis */
	normalize (NormalAxes[i]->direction);
	for (k=0; k<DIMENSION; k++) norm[k] = NormalAxes[i]->direction[k];
	i_one = -1;
	n_one = n_zero = 0;
	for (k=0; k<DIMENSION; k++)
	{
	  if (fabs(fabs(norm[k])-1) <= ToleranceFinal)
	  {
	    i_one = k;
	    n_one ++;
	  }
	  if (fabs(norm[k]) <= ToleranceFinal) n_zero ++;
	}
	if (n_one!=1 || n_zero!=2) continue;		/* not aligned axis */
	if (caxes[i_one])          continue; /*this direction already counted*/
	if (i_one == 2)		   continue;	/* not perpendicular axis */
	naxes ++;
	caxes[i_one] = 1;					/* count it */
      }
      if (naxes == 2) return;
      if (naxes == 1)			/* align single order 2 axis along X */
      {
	i_one = 0;
	for (k=0; k<DIMENSION; k++) if (caxes[k]) i_one = k;   /* find coord */
	swap_x_coords (i_one);
	return;							/* Dnd */
      }
      /* Align parallel symmetry planes (if any) */
      nplanes = 0;
      for (k=0; k<DIMENSION; k++) cplanes[k] = 0;
      for (i=0; i<PlanesCount; i++)
      {
	normalize (Planes[i]->normal);
	for (k=0; k<DIMENSION; k++) norm[k] = Planes[i]->normal[k];
	i_one = -1;
	n_one = n_zero = 0;
	for (k=0; k<DIMENSION; k++)
	{
	  if (fabs(fabs(norm[k])-1) <= ToleranceFinal)
	  {
	    i_one = k;
	    n_one ++;
	  }
	  if (fabs(norm[k]) <= ToleranceFinal) n_zero ++;
	}
	if (n_one!=1 || n_zero!=2) continue;		/* not aligned plane */
	if (cplanes[i_one])        continue; /*this direction already counted*/
	if (i_one == 2)		   continue;	/* perpendicular plane */
	nplanes ++;
	cplanes[i_one] = 1;					/* count it */
      }
      if (nplanes == 2) return;
      if (nplanes == 1)		/* align single plane perpendicular to Y */
      {
	i_one = 1;
	for (k=0; k<DIMENSION; k++) if (cplanes[k]) i_one = k; /* find coord */
	swap_y_coords (i_one);
	return;							/* Cnv */
      }
      return;
    }
  }

  /* Align symmetry planes (if any) */
  nplanes = 0;
  for (k=0; k<DIMENSION; k++) cplanes[k] = 0;
  for (i=0; i<PlanesCount; i++)
  {
    normalize (Planes[i]->normal);
    for (k=0; k<DIMENSION; k++) norm[k] = Planes[i]->normal[k];
    i_one = -1;
    n_one = n_zero = 0;
    for (k=0; k<DIMENSION; k++)
    {
      if (fabs(fabs(norm[k])-1) <= ToleranceFinal)
      {
	i_one = k;
	n_one ++;
      }
      if (fabs(norm[k]) <= ToleranceFinal) n_zero ++;
    }
    if (n_one!=1 || n_zero!=2) continue;		/* not abelian plane */
    if (cplanes[i_one])        continue;   /* this direction already counted */
    nplanes ++;
    cplanes[i_one] = 1;						/* count it */
  }
  if (nplanes == 3) return;	/* D2h, all directions used, nothing to swap */
  if (nplanes == 2)				/* align two planes along Z */
  {
    i_zero = 2;
    for (k=0; k<DIMENSION; k++) if (! cplanes[k]) i_zero = k;  /* find coord */
    swap_coords (i_zero);
    return;							/* C2v */
  }
  if (nplanes == 1)		/* align single plane perpendicular to Z */
  {
    i_one = 2;
    for (k=0; k<DIMENSION; k++) if (cplanes[k]) i_one = k;     /* find coord */
    swap_coords (i_one);
    return;							/* Cs or C2h */
  }

  /* Align order 2 normal axes (if any) */
  naxes = 0;
  for (k=0; k<DIMENSION; k++) caxes[k] = 0;
  for (i=0; i<NormalAxesCount; i++)
  {
    if (NormalAxes[i]->order != 2) continue;	/* not an order 2 axis */
    normalize (NormalAxes[i]->direction);
    for (k=0; k<DIMENSION; k++) norm[k] = NormalAxes[i]->direction[k];
    i_one = -1;
    n_one = n_zero = 0;
    for (k=0; k<DIMENSION; k++)
    {
      if (fabs(fabs(norm[k])-1) <= ToleranceFinal)
      {
	i_one = k;
	n_one ++;
      }
      if (fabs(norm[k]) <= ToleranceFinal) n_zero ++;
    }
    if (n_one!=1 || n_zero!=2) continue;		/* not abelian axis */
    if (caxes[i_one])          continue;   /* this direction already counted */
    naxes ++;
    caxes[i_one] = 1;						/* count it */
  }
  if (naxes == 3) return;	/* D2, all directions used, nothing to swap */
  if (naxes == 2)		/* align two order 2 axes perpendicular to Z */
  {
    i_zero = 2;
    for (k=0; k<DIMENSION; k++) if (! caxes[k]) i_zero = k;    /* find coord */
    swap_coords (i_zero);
    return;		/* no planes, two order 2 axes -- is it possible ? */
  }
  if (naxes == 1)			/* align single order 2 axis along Z */
  {
    i_one = 2;
    for (k=0; k<DIMENSION; k++) if (caxes[k]) i_one = k;       /* find coord */
    swap_coords (i_one);
    return;							/* C2 */
  }

  /* Align normal axes of any order */
  naxes = 0;
  for (k=0; k<DIMENSION; k++) caxes[k] = 0;
  for (i=0; i<NormalAxesCount; i++)
  {
    normalize (NormalAxes[i]->direction);
    for (k=0; k<DIMENSION; k++) norm[k] = NormalAxes[i]->direction[k];
    i_one = -1;
    n_one = n_zero = 0;
    for (k=0; k<DIMENSION; k++)
    {
      if (fabs(fabs(norm[k])-1) <= ToleranceFinal)
      {
	i_one = k;
	n_one ++;
      }
      if (fabs(norm[k]) <= ToleranceFinal) n_zero ++;
    }
    if (n_one!=1 || n_zero!=2) continue;	/* not along coordinate axis */
    if (caxes[i_one])          continue;   /* this direction already counted */
    naxes ++;
    caxes[i_one] = 1;						/* count it */
  }
  if (naxes == 3) return;	/* all directions used -- nothing to swap ? */
  if (naxes == 2)			/* align two axes perpendicular to Z */
  {
    i_zero = 2;
    for (k=0; k<DIMENSION; k++) if (! caxes[k]) i_zero = k;    /* find coord */
    swap_coords (i_zero);
    return;	/* no planes, no order 2 axis, two axes -- is it possible ? */
  }
  if (naxes == 1)				/* align single axis along Z */
  {
    i_one = 2;
    for (k=0; k<DIMENSION; k++) if (caxes[k]) i_one = k;       /* find coord */
    swap_coords (i_one);
    return;						/* C3, C5, C7, S6 */
  }

  /* Align improper axes (no planes, no normal axes -- is it possible ?) */
  naxes = 0;
  for (k=0; k<DIMENSION; k++) caxes[k] = 0;
  for (i=0; i<ImproperAxesCount; i++)
  {
    normalize (ImproperAxes[i]->direction);
    for (k=0; k<DIMENSION; k++) norm[k] = ImproperAxes[i]->direction[k];
    i_one = -1;
    n_one = n_zero = 0;
    for (k=0; k<DIMENSION; k++)
    {
      if (fabs(fabs(norm[k])-1) <= ToleranceFinal)
      {
	i_one = k;
	n_one ++;
      }
      if (fabs(norm[k]) <= ToleranceFinal) n_zero ++;
    }
    if (n_one!=1 || n_zero!=2) continue;	/* not along coordinate axis */
    if (caxes[i_one])          continue;   /* this direction already counted */
    naxes ++;
    caxes[i_one] = 1;						/* count it */
  }
  if (naxes == 3) return;	/* all directions used -- nothing to swap */
  if (naxes == 2)	/* align two improper axes perpendicular to Z */
  {
    i_zero = 2;
    for (k=0; k<DIMENSION; k++) if (! caxes[k]) i_zero = k;    /* find coord */
    swap_coords (i_zero);
    return;
  }
  if (naxes == 1)		/* align single improper axis along Z */
  {
    i_one = 2;
    for (k=0; k<DIMENSION; k++) if (caxes[k]) i_one = k;       /* find coord */
    swap_coords (i_one);
    return;
  }
}

/*************************************************************************
 * Symmetrize atoms around either all or only abelian symmetry planes
 * to eliminate distortions due to roundoff errors.
 * For each atom in sequence another atom is searched for which is closest
 * to its mirror. If closest to the mirror is the original atom itself
 * it is averaged with its own mirror. If some another atom is closest
 * each atom of the pair is averaged with the mirror of the other.
 * The mirror atom is marked in order not to take it later as a mirror for
 * some third atom.
 *************************************************************************/
static void symmetrize_planes (int check_abel)
{
  int		i, j, k, n, i_min, n_one, n_zero;
  double	dmin, dist;
  ATOM1		atom;
  
  for (i=0; i<PlanesCount; i++)				/* for all planes */
  {
    if ((SymmPlane || SymmAxis || SymmImproper || SymmCenter) &&
	i != SymmPlane-1)
      continue;		/* dont symmetrize around this particular plane */

    if (check_abel)				/* abelian planes only ? */
    {
      if (fabs(Planes[i]->distance) > ToleranceFinal) continue;	/*not abelian*/
      n_one = n_zero = 0;
      for (k=0; k<DIMENSION; k++)
      {
	if (fabs(fabs(Planes[i]->normal[k])-1) <= ToleranceFinal) n_one  ++;
	if (fabs(Planes[i]->normal[k])         <= ToleranceFinal) n_zero ++;
      }
      if (n_one!=1 || n_zero!=2) continue;		/* not abelian */
    }

    for (j=0; j<AtomsCount; j++) Atoms[j].type = abs(Atoms[j].type); /*unmark*/

    for (j=0; j<AtomsCount; j++)			/* for all atoms */
    {
      if (Atoms[j].type < 0) continue;	/* already marked as some mirror */
      mirror_atom (Planes[i], Atoms+j, &atom);
      dmin = 0;
      for (k=0; k<DIMENSION; k++)
	dmin += (Atoms[j].x[k]-atom.x[k])*(Atoms[j].x[k]-atom.x[k]);
      i_min = j;
      for (n=j+1; n<AtomsCount; n++)	/* search the closest unmarked atom */
      {
	if (Atoms[n].type < 0)          continue;	/* already marked */
	if (Atoms[n].type != atom.type) continue;	/* wrong atom number */
	dist = 0;
	for (k=0; k<DIMENSION; k++)
	  dist += (Atoms[n].x[k]-atom.x[k])*(Atoms[n].x[k]-atom.x[k]);
	if (dist < dmin)		/* this one is closer, memorize it */
	{
	  dmin  = dist;
	  i_min = n;
	}
      }
      for (k=0; k<DIMENSION; k++)
	Atoms[i_min].x[k] = (Atoms[i_min].x[k]+atom.x[k])/2;	/* average */
      if (i_min == j) continue;			/* average with itself */
      mirror_atom (Planes[i], Atoms+i_min, Atoms+j);	/* reflect average */
      Atoms[i_min].type = -Atoms[j].type;			/* and mark */
    }							/* for all atoms */
  }							/* for all planes */
}

/*************************************************************************
 * Symmetrize atoms around either all or only abelian normal axes
 * to eliminate distortions due to roundoff errors.
 * For each atom in sequence other atoms are searched for which are closest
 * to all its images. If closest to the images is the original atom itself
 * it is averaged with all its images. If some other atoms are closest
 * each atom of the set is averaged with the images of all others.
 * The image atoms are marked in order not to take them later as images for
 * some third atom.
 * Analogous to symmetrize_planes() but complicated by the fact that there
 * may be more or less than one rotation images depending on the axis order.
 *************************************************************************/
static void symmetrize_axes (int check_abel)
{
  int		i, j, k, m, n, n_sym, i_min, n_one, n_zero, order;
  double	dmin, dist;
  ATOM1		atom, atom2, mean;
  
  for (i=0; i<NormalAxesCount; i++)		/* for all normal axes */
  {
    if ((SymmPlane || SymmAxis || SymmImproper || SymmCenter) &&
	i != SymmAxis-1)
      continue;		/* dont symmetrize around this particular axis */

    if (check_abel)			/* abelian normal axes only ? */
    {					/* only order 2 or inf allowed */
      if (NormalAxes[i]->order && NormalAxes[i]->order != 2) continue;
      if (fabs(NormalAxes[i]->distance) > ToleranceFinal)    continue;/* bad */
      n_one = n_zero = 0;
      for (k=0; k<DIMENSION; k++)
      {
	if (fabs(fabs(NormalAxes[i]->direction[k])-1) <= ToleranceFinal)
	  n_one ++;
	if (fabs(NormalAxes[i]->direction[k]) <= ToleranceFinal) n_zero ++;
      }
      if (n_one!=1 || n_zero!=2) continue;		/* not abelian */
    }

    order = NormalAxes[i]->order;		/* memorize original order */
    if (! NormalAxes[i]->order) NormalAxes[i]->order = 2;/* make inf order 2 */
    for (j=0; j<AtomsCount; j++) Atoms[j].type = abs(Atoms[j].type); /*unmark*/

    for (j=0; j<AtomsCount; j++)			/* for all atoms */
    {
      if (Atoms[j].type < 0) continue;	/* already marked as some image */
      rotate_atom (NormalAxes[i], Atoms+j, &atom);
      if (! order)	/* put atom on the axis if it has infinity order */
      {
	for (k=0; k<DIMENSION; k++)
	  Atoms[j].x[k] = (Atoms[j].x[k]+atom.x[k])/2;
	continue;	/* nothing to do except averaging with itself */
      }
      n_sym = 1;
      dmin  = 0;
      for (k=0; k<DIMENSION; k++)
	dmin += (Atoms[j].x[k]-atom.x[k])*(Atoms[j].x[k]-atom.x[k]);
      i_min = j;
      for (n=j+1; n<AtomsCount; n++)	/* search the closest unmarked atom */
      {
	if (Atoms[n].type < 0)          continue;	/* already marked */
	if (Atoms[n].type != atom.type) continue;	/* wrong atom number */
	dist = 0;
	for (k=0; k<DIMENSION; k++)
	  dist += (Atoms[n].x[k]-atom.x[k])*(Atoms[n].x[k]-atom.x[k]);
	if (dist < dmin)		/* this one is closer, memorize it */
	{
	  dmin  = dist;
	  i_min = n;
	}
      }
      if (i_min == j)				/* average with itself */
      {
	for (k=0; k<DIMENSION; k++) mean.x[k] = Atoms[j].x[k]+atom.x[k];
	for (m=2; m<order; m++)		/* take more images into account */
	{
	  rotate_atom (NormalAxes[i], &atom, &atom2);
	  for (k=0; k<DIMENSION; k++)
	  {
	    atom.x[k] = atom2.x[k];
	    mean.x[k] += atom.x[k];
	  }
	}
	for (k=0; k<DIMENSION; k++) Atoms[j].x[k] = mean.x[k]/order;
	continue;					/* and nothing else */
      }

      Atoms[i_min].type = 10000+n_sym++;		/* count this image */
      for (m=2; m<order; m++)				/* count more images */
      {
	rotate_atom (NormalAxes[i], &atom, &atom2);
	for (k=0; k<DIMENSION; k++) atom.x[k] = atom2.x[k];
	dmin = 0;
	for (k=0; k<DIMENSION; k++)
	  dmin += (Atoms[j].x[k]-atom.x[k])*(Atoms[j].x[k]-atom.x[k]);
	i_min = j;
	for (n=j+1; n<AtomsCount; n++)	/* search the closest uncounted atom */
	{		/* reject wrong number, marked or already counted */
	  if (Atoms[n].type < 0)                                    continue;
	  if (Atoms[n].type != atom.type && Atoms[n].type <= 10000) continue;
	  dist = 0;
	  for (k=0; k<DIMENSION; k++)
	    dist += (Atoms[n].x[k]-atom.x[k])*(Atoms[n].x[k]-atom.x[k]);
	  if (dist < dmin)		/* this one is closer, memorize it */
	  {
	    dmin  = dist;
	    i_min = n;
	  }
	}
	if (i_min == j || Atoms[i_min].type > 10000) continue;/*nothing found*/
	Atoms[i_min].type = 10000+n_sym++;    /* count this additional image */
      }								/* for order */
      if (n_sym != order)      /* number of images does not match axis order */
      {
	for (n=j+1; n<AtomsCount; n++)		/* uncount counted atoms */
	  if (Atoms[n].type > 10000) Atoms[n].type = Atoms[j].type;
	rotate_atom (NormalAxes[i], Atoms+j, &atom);
	for (k=0; k<DIMENSION; k++) mean.x[k] = Atoms[j].x[k]+atom.x[k];
	for (m=2; m<order; m++)			/* do average with itself */
	{
	  rotate_atom (NormalAxes[i], &atom, &atom2);
	  for (k=0; k<DIMENSION; k++)
	  {
	    atom.x[k] = atom2.x[k];
	    mean.x[k] += atom.x[k];
	  }
	}
	for (k=0; k<DIMENSION; k++) Atoms[j].x[k] = mean.x[k]/order;
	continue;					/* nothing else */
      }

      for (k=0; k<DIMENSION; k++) mean.x[k] = Atoms[j].x[k];	/* images OK */
      for (n=j+1; n<AtomsCount; n++)	/* average with all counted images */
      {
	n_sym = Atoms[n].type-10000;	/* extract number of rotations */
	if (n_sym<=0 || n_sym>=order) continue;	/* not a counted image */
	rotate_atom (NormalAxes[i], Atoms+n, &atom);
	for (m=n_sym+1; m<order; m++)	/* accumulate proper rotation angle */
	{
	  rotate_atom (NormalAxes[i], &atom, &atom2);
	  for (k=0; k<DIMENSION; k++) atom.x[k] = atom2.x[k];
	}
	for (k=0; k<DIMENSION; k++) mean.x[k] += atom.x[k];
      }
      for (k=0; k<DIMENSION; k++) Atoms[j].x[k] = mean.x[k]/order;
      for (n=j+1; n<AtomsCount; n++)	/* make proper number of images */
      {
	n_sym = Atoms[n].type-10000;	/* extract number of rotations */
	if (n_sym<=0 || n_sym>=order) continue;	/* not a counted image */
	rotate_atom (NormalAxes[i], Atoms+j, &atom);
	for (m=1; m<n_sym; m++)		/* accumulate proper rotation angle */
	{
	  rotate_atom (NormalAxes[i], &atom, &atom2);
	  for (k=0; k<DIMENSION; k++) atom.x[k] = atom2.x[k];
	}
	for (k=0; k<DIMENSION; k++) Atoms[n].x[k] = atom.x[k];
	Atoms[n].type = -Atoms[j].type;				/* and mark */
      }
    }							/* for all atoms */
    NormalAxes[i]->order = order;		/* restore original order */
  }						/* for all normal axes */
}

/*************************************************************************
 * Symmetrize atoms around either all or only abelian improper axes
 * to eliminate distortions due to roundoff errors.
 * For each atom in sequence other atoms are searched for which are closest
 * to all its images. If closest to the images is the original atom itself
 * it is averaged with all its images. If some other atoms are closest
 * each atom of the set is averaged with the images of all others.
 * The image atoms are marked in order not to take them later as images for
 * some third atom.
 * Analogous to symmetrize_axes() but complicated by the fact
 * that the reflection plane has to be taken into account
 * even for the atoms laying on the improper axis.
 *************************************************************************/
static void symmetrize_improper_axes (int check_abel)
{
  int		i, j, k, m, n, n_sym, i_min, n_one, n_zero, order;
  double	dmin, dist;
  ATOM1		atom, atom2, mean;
  
  for (i=0; i<ImproperAxesCount; i++)		/* for all improper axes */
  {
    if ((SymmPlane || SymmAxis || SymmImproper || SymmCenter) &&
	i != SymmImproper-1)
      continue;	/* dont symmetrize around this particular improper axis */

    if (check_abel)			/* abelian improper axes only ? */
    {						/* only order 2 allowed */
      if (ImproperAxes[i]->order != 2)                      continue;
      if (fabs(ImproperAxes[i]->distance) > ToleranceFinal) continue; /* bad */
      n_one = n_zero = 0;
      for (k=0; k<DIMENSION; k++)
      {
	if (fabs(fabs(ImproperAxes[i]->direction[k])-1) <= ToleranceFinal)
	  n_one ++;
	if (fabs(ImproperAxes[i]->direction[k]) <= ToleranceFinal) n_zero ++;
      }
      if (n_one!=1 || n_zero!=2) continue;		/* not abelian */
    }

    for (j=0; j<AtomsCount; j++) Atoms[j].type = abs(Atoms[j].type); /*unmark*/

    for (j=0; j<AtomsCount; j++)			/* for all atoms */
    {
      if (Atoms[j].type < 0) continue;	/* already marked as some image */
      rotate_reflect_atom (ImproperAxes[i], Atoms+j, &atom);
      n_sym = 1;
      dmin  = 0;
      for (k=0; k<DIMENSION; k++)
	dmin += (Atoms[j].x[k]-atom.x[k])*(Atoms[j].x[k]-atom.x[k]);
      i_min = j;
      for (n=j+1; n<AtomsCount; n++)	/* search the closest unmarked atom */
      {
	if (Atoms[n].type < 0)          continue;	/* already marked */
	if (Atoms[n].type != atom.type) continue;	/* wrong atom number */
	dist = 0;
	for (k=0; k<DIMENSION; k++)
	  dist += (Atoms[n].x[k]-atom.x[k])*(Atoms[n].x[k]-atom.x[k]);
	if (dist < dmin)		/* this one is closer, memorize it */
	{
	  dmin  = dist;
	  i_min = n;
	}
      }
      if (i_min == j)				/* average with itself */
      {
	rotate_atom (ImproperAxes[i], Atoms+j, &atom);
	for (k=0; k<DIMENSION; k++) mean.x[k] = Atoms[j].x[k]+atom.x[k];
	for (m=2; m<ImproperAxes[i]->order; m++)      /* account more images */
	{
	  rotate_atom (ImproperAxes[i], &atom, &atom2);
	  for (k=0; k<DIMENSION; k++)
	  {
	    atom.x[k] = atom2.x[k];
	    mean.x[k] += atom.x[k];
	  }
	}
	for (k=0; k<DIMENSION; k++)	/* average without reflection */
	  Atoms[j].x[k] = mean.x[k]/ImproperAxes[i]->order;
	rotate_reflect_atom (ImproperAxes[i], Atoms+j, &atom);
	for (k=0; k<DIMENSION; k++)		/* and with reflection */
	  Atoms[j].x[k] = (Atoms[j].x[k]+atom.x[k])/2;
	continue;					/* and nothing else */
      }

      Atoms[i_min].type = 10000+n_sym++;		/* count this image */
      order = ImproperAxes[i]->order;
      if (order & 1) order <<= 1;
      for (m=2; m<order; m++)				/* count more images */
      {
	rotate_reflect_atom (ImproperAxes[i], &atom, &atom2);
	for (k=0; k<DIMENSION; k++) atom.x[k] = atom2.x[k];
	dmin = 0;
	for (k=0; k<DIMENSION; k++)
	  dmin += (Atoms[j].x[k]-atom.x[k])*(Atoms[j].x[k]-atom.x[k]);
	i_min = j;
	for (n=j+1; n<AtomsCount; n++)	/* search the closest uncounted atom */
	{		/* reject wrong number, marked or already counted */
	  if (Atoms[n].type < 0)                                    continue;
	  if (Atoms[n].type != atom.type && Atoms[n].type <= 10000) continue;
	  dist = 0;
	  for (k=0; k<DIMENSION; k++)
	    dist += (Atoms[n].x[k]-atom.x[k])*(Atoms[n].x[k]-atom.x[k]);
	  if (dist < dmin)		/* this one is closer, memorize it */
	  {
	    dmin  = dist;
	    i_min = n;
	  }
	}
	if (i_min == j || Atoms[i_min].type > 10000) continue;/*nothing found*/
	Atoms[i_min].type = 10000+n_sym++;    /* count this additional image */
      }								/* for order */

      if (n_sym == order)	/* complete number of images available */
      {
	for (k=0; k<DIMENSION; k++) mean.x[k] = Atoms[j].x[k];	/* images OK */
	for (n=j+1; n<AtomsCount; n++)	/* average with all counted images */
	{
	  n_sym = Atoms[n].type-10000;	/* extract number of transformations */
	  if (n_sym<=0 || n_sym>=order) continue;	/* not counted */
	  rotate_reflect_atom (ImproperAxes[i], Atoms+n, &atom);
	  for (m=n_sym+1; m<order; m++)	/* accumulate transformations */
	  {
	    rotate_reflect_atom (ImproperAxes[i], &atom, &atom2);
	    for (k=0; k<DIMENSION; k++) atom.x[k] = atom2.x[k];
	  }
	  for (k=0; k<DIMENSION; k++) mean.x[k] += atom.x[k];
	}
	for (k=0; k<DIMENSION; k++) Atoms[j].x[k] = mean.x[k]/order;
	for (n=j+1; n<AtomsCount; n++)	/* make proper number of images */
	{
	  n_sym = Atoms[n].type-10000;	/* extract number of transformations */
	  if (n_sym<=0 || n_sym>=order) continue;	/* not counted */
	  rotate_reflect_atom (ImproperAxes[i], Atoms+j, &atom);
	  for (m=1; m<n_sym; m++)  /* accumulate proper number of transforms */
	  {
	    rotate_reflect_atom (ImproperAxes[i], &atom, &atom2);
	    for (k=0; k<DIMENSION; k++) atom.x[k] = atom2.x[k];
	  }
	  for (k=0; k<DIMENSION; k++) Atoms[n].x[k] = atom.x[k];
	  Atoms[n].type = -Atoms[j].type;			/* and mark */
	}
	continue;					/* nothing else */
      }				/* if complete number of images available */

      if (n_sym == ImproperAxes[i]->order)  /* atoms lie on reflection plane */
      {
	for (k=0; k<DIMENSION; k++) mean.x[k] = Atoms[j].x[k];	/* images OK */
	for (n=j+1; n<AtomsCount; n++)	/* average with all counted images */
	{
	  n_sym = Atoms[n].type-10000;	/* extract number of transformations */
	  if (n_sym<=0 || n_sym>=ImproperAxes[i]->order) continue;/*uncounted*/
	  rotate_atom (ImproperAxes[i], Atoms+n, &atom);
	  for (m=n_sym+1; m<ImproperAxes[i]->order; m++)       /* accumulate */
	  {
	    rotate_atom (ImproperAxes[i], &atom, &atom2);
	    for (k=0; k<DIMENSION; k++) atom.x[k] = atom2.x[k];
	  }
	  for (k=0; k<DIMENSION; k++) mean.x[k] += atom.x[k];
	}
	for (k=0; k<DIMENSION; k++)	/* average with rotated images */
	  Atoms[j].x[k] = mean.x[k]/ImproperAxes[i]->order;
	rotate_reflect_atom (ImproperAxes[i], Atoms+j, &atom);
	for (m=1; m<ImproperAxes[i]->order; m++) /*accumulate transformations*/
	{
	  rotate_reflect_atom (ImproperAxes[i], &atom, &atom2);
	  for (k=0; k<DIMENSION; k++) atom.x[k] = atom2.x[k];
	}
	for (k=0; k<DIMENSION; k++)	/* average with reflected itself */
	  Atoms[j].x[k] = (Atoms[j].x[k]+atom.x[k])/2;
	for (n=j+1; n<AtomsCount; n++)	/* make proper number of images */
	{
	  n_sym = Atoms[n].type-10000;	/* extract number of transformations */
	  if (n_sym<=0 || n_sym>=ImproperAxes[i]->order) continue;/*uncounted*/
	  rotate_atom (ImproperAxes[i], Atoms+j, &atom);
	  for (m=1; m<n_sym; m++)  /* accumulate proper number of transforms */
	  {
	    rotate_atom (ImproperAxes[i], &atom, &atom2);
	    for (k=0; k<DIMENSION; k++) atom.x[k] = atom2.x[k];
	  }
	  for (k=0; k<DIMENSION; k++) Atoms[n].x[k] = atom.x[k];
	  Atoms[n].type = -Atoms[j].type;			/* and mark */
	}
	continue;					/* nothing else */
      }					/* if atoms lie on reflection plane */

      if (n_sym == 2)			/* atoms lie on rotation axis */
      {			/* the reflected image has to be averaged separately */
	rotate_atom (ImproperAxes[i], Atoms+j, &atom);
	for (k=0; k<DIMENSION; k++) mean.x[k] = Atoms[j].x[k]+atom.x[k];
	for (m=2; m<ImproperAxes[i]->order; m++)     /* accumulate rotations */
	{
	  rotate_atom (ImproperAxes[i], &atom, &atom2);
	  for (k=0; k<DIMENSION; k++)
	  {
	    atom.x[k] = atom2.x[k];
	    mean.x[k] += atom.x[k];
	  }
	}
	for (k=0; k<DIMENSION; k++)	/* average with rotated itself */
	  Atoms[j].x[k] = mean.x[k]/ImproperAxes[i]->order;
	for (n=j+1; n<AtomsCount; n++) if (Atoms[n].type == 10001) break;
	rotate_atom (ImproperAxes[i], Atoms+n, &atom);
	for (k=0; k<DIMENSION; k++) mean.x[k] = Atoms[n].x[k]+atom.x[k];
	for (m=2; m<ImproperAxes[i]->order; m++)     /* accumulate rotations */
	{
	  rotate_atom (ImproperAxes[i], &atom, &atom2);
	  for (k=0; k<DIMENSION; k++)
	  {
	    atom.x[k] = atom2.x[k];
	    mean.x[k] += atom.x[k];
	  }
	}
	for (k=0; k<DIMENSION; k++)    /* average mirror with rotated itself */
	  Atoms[n].x[k] = mean.x[k]/ImproperAxes[i]->order;
	rotate_reflect_atom (ImproperAxes[i], Atoms+n, &atom);
	for (k=0; k<DIMENSION; k++)	/* average with reflected mirror */
	  Atoms[j].x[k] = (Atoms[j].x[k]+atom.x[k])/2;
	rotate_reflect_atom (ImproperAxes[i], Atoms+j, Atoms+n);  /* reflect */
	Atoms[n].type = -Atoms[j].type;				/* and mark */
	continue;					/* nothing else */
      }					/* if atoms lie on rotation axis */
    }							/* for all atoms */
  }						/* for all improper axes */
}

/*************************************************************************
 * Symmetrize atoms around either all or only abelian inversion centers
 * to eliminate distortions due to roundoff errors.
 * For each atom in sequence another atom is searched for which is closest
 * to its inversion. If closest to the inverted atom is the original atom
 * itself it is averaged with its own inversion. If some another atom is
 * closest each atom of the pair is averaged with the inversion of the other.
 * The inverted atom is marked in order not to take it later as an inversion
 * for some third atom.
 * Similar to symmetrize_planes().
 *************************************************************************/
static void symmetrize_inversion_centers (int check_abel)
{
  int		i, j, k, n, i_min;
  double	dmin, dist;
  ATOM1		atom;
  
  for (i=0; i<InversionCentersCount; i++)		/* for all centers */
  {
    if ((SymmPlane || SymmAxis || SymmImproper || SymmCenter) &&
	i != SymmCenter-1)
      continue;	/* dont symmetrize around this particular inversion center */

    if (check_abel && fabs(InversionCenters[i]->distance) > ToleranceFinal)
      continue;						/* not abelian */

    for (j=0; j<AtomsCount; j++) Atoms[j].type = abs(Atoms[j].type); /*unmark*/

    for (j=0; j<AtomsCount; j++)			/* for all atoms */
    {
      if (Atoms[j].type < 0) continue;	/* already marked as some inversion */
      invert_atom (InversionCenters[i], Atoms+j, &atom);
      dmin = 0;
      for (k=0; k<DIMENSION; k++)
	dmin += (Atoms[j].x[k]-atom.x[k])*(Atoms[j].x[k]-atom.x[k]);
      i_min = j;
      for (n=j+1; n<AtomsCount; n++)	/* search the closest unmarked atom */
      {
	if (Atoms[n].type < 0)          continue;	/* already marked */
	if (Atoms[n].type != atom.type) continue;	/* wrong atom number */
	dist = 0;
	for (k=0; k<DIMENSION; k++)
	  dist += (Atoms[n].x[k]-atom.x[k])*(Atoms[n].x[k]-atom.x[k]);
	if (dist < dmin)		/* this one is closer, memorize it */
	{
	  dmin  = dist;
	  i_min = n;
	}
      }
      for (k=0; k<DIMENSION; k++)
	Atoms[i_min].x[k] = (Atoms[i_min].x[k]+atom.x[k])/2;	/* average */
      if (i_min == j) continue;			/* average with itself */
      invert_atom (InversionCenters[i], Atoms+i_min, Atoms+j);	/* invert */
      Atoms[i_min].type = -Atoms[j].type;			/* and mark */
    }							/* for all atoms */
  }							/* for all centers */
}

/*************************************************************************
 * Mark atoms symmetric around either all or only abelian planes
 * to find out symmetry unique atoms.
 * For each atom in sequence another atom is searched for which is
 * sufficiently close to its mirror. If it is close enough, it is marked.
 * There is some similarity with symmetrize_planes()
 * but the criterium of closeness is different.
 *************************************************************************/
static void mark_symmetric_planes (int check_abel)
{
  int	i, j, k, n, n_one, n_zero;
  ATOM1	atom;
  
  for (i=0; i<PlanesCount; i++)				/* for all planes */
  {
    if (check_abel)				/* abelian planes only ? */
    {
      if (fabs(Planes[i]->distance) > ToleranceFinal) continue;	/*not abelian*/
      n_one = n_zero = 0;
      for (k=0; k<DIMENSION; k++)
      {
	if (fabs(fabs(Planes[i]->normal[k])-1) <= ToleranceFinal) n_one  ++;
	if (fabs(Planes[i]->normal[k])         <= ToleranceFinal) n_zero ++;
      }
      if (n_one!=1 || n_zero!=2) continue;		/* not abelian */
    }

    for (j=0; j<AtomsCount; j++)			/* for all atoms */
    {
      if (Atoms[j].type < 0) continue;	/* already marked as non unique */
      mirror_atom (Planes[i], Atoms+j, &atom);
      n_zero = 0;
      for (k=0; k<DIMENSION; k++)
	if (fabs(Atoms[j].x[k]-atom.x[k]) <= ToleranceFinal) n_zero ++;
      if (n_zero == 3) continue;      /* the original itself is close enough */
      for (n=j+1; n<AtomsCount; n++)   /* search another close unmarked atom */
      {
	if (Atoms[n].type < 0) continue;     /* already marked as non unique */
	n_zero = 0;
	for (k=0; k<DIMENSION; k++)
	  if (fabs(Atoms[n].x[k]-atom.x[k]) <= ToleranceFinal) n_zero ++;
	if (Atoms[n].type == atom.type) n_zero ++;	/* check atom number */
	if (n_zero == 4)			/* this one is close enough */
	{
	  Atoms[n].type = -Atoms[j].type;			/* mark it */
	  break;
	}
      }
    }							/* for all atoms */
  }							/* for all planes */
}

/*************************************************************************
 * Mark atoms symmetric around either all or only abelian normal axes
 * to find out symmetry unique atoms.
 * For each atom in sequence other atoms are searched for which are
 * sufficiently close to all its rotation images.
 * If they are close enough, they are marked.
 * There is some similarity with symmetrize_axes()
 * but the criterium of closeness is different.
 *************************************************************************/
static void mark_symmetric_axes (int check_abel)
{
  int	i, j, k, m, n, n_one, n_zero;
  ATOM1	atom, atom2;
  
  for (i=0; i<NormalAxesCount; i++)		/* for all normal axes */
  {
    if (! NormalAxes[i]->order) continue;    /* nothing to mark for infinity */
    if (check_abel)			/* abelian normal axes only ? */
    {
      if (NormalAxes[i]->order != 2) continue;	/* only order 2 allowed */
      if (fabs(NormalAxes[i]->distance) > ToleranceFinal) continue;   /* bad */
      n_one = n_zero = 0;
      for (k=0; k<DIMENSION; k++)
      {
	if (fabs(fabs(NormalAxes[i]->direction[k])-1) <= ToleranceFinal)
	  n_one ++;
	if (fabs(NormalAxes[i]->direction[k]) <= ToleranceFinal) n_zero ++;
      }
      if (n_one!=1 || n_zero!=2) continue;		/* not abelian */
    }

    for (j=0; j<AtomsCount; j++)			/* for all atoms */
    {
      if (Atoms[j].type < 0) continue;	/* already marked as non unique */
      rotate_atom (NormalAxes[i], Atoms+j, &atom);
      n_zero = 0;
      for (k=0; k<DIMENSION; k++)
	if (fabs(Atoms[j].x[k]-atom.x[k]) <= ToleranceFinal) n_zero ++;
      if (n_zero == 3) continue;      /* the original itself is close enough */
      for (m=1; m<NormalAxes[i]->order; m++) /*accumulate all rotation angles*/
      {
	for (n=j+1; n<AtomsCount; n++) /* search another close unmarked atom */
	{
	  if (Atoms[n].type < 0) continue;   /* already marked as non unique */
	  n_zero = 0;
	  for (k=0; k<DIMENSION; k++)
	    if (fabs(Atoms[n].x[k]-atom.x[k]) <= ToleranceFinal) n_zero ++;
	  if (Atoms[n].type == atom.type) n_zero ++;	/* check atom number */
	  if (n_zero == 4)			/* this one is close enough */
	  {
	    Atoms[n].type = -Atoms[j].type;			/* mark it */
	    break;
	  }
	}
	rotate_atom (NormalAxes[i], &atom, &atom2);
	for (k=0; k<DIMENSION; k++) atom.x[k] = atom2.x[k];
      }
    }							/* for all atoms */
  }						/* for all normal axes */
}
  
/*************************************************************************
 * Mark atoms symmetric around principal axes
 * to find out symmetry unique atoms.
 * For each atom in sequence other atoms are searched for which are
 * sufficiently close to all its rotation images.
 * If they are close enough, they are marked.
 * There is some similarity with mark_symmetric_axes()
 * but the unique atoms lying in the YZ or XZ planes are preferred.
 *************************************************************************/
static void mark_principal_axes (void)
{
  int		i, j, k, m, n, n_one, n_zero, maxis, i_one;
  ATOM1		atom, atom2;
  double	norm[DIMENSION];

  /* Find principal axis order */
  maxis = 0;
  for (i=0; i<NormalAxesCount; i++)
  {
    if (! NormalAxes[i]->order) continue;		/* an infinity axis */
    normalize (NormalAxes[i]->direction);
    for (k=0; k<DIMENSION; k++) norm[k] = NormalAxes[i]->direction[k];
    n_one = n_zero = 0;
    for (k=0; k<DIMENSION; k++)
    {
      if (fabs(fabs(norm[k])-1) <= ToleranceFinal) n_one  ++;
      if (fabs(norm[k])         <= ToleranceFinal) n_zero ++;
    }
    if (n_one!=1 || n_zero!=2) continue;		/* not aligned axis */
    if (maxis < NormalAxes[i]->order) maxis = NormalAxes[i]->order;
  }
  if (maxis <= 2) return;			/* uninteresting axis order */

  /* Operate with all principal axes */
  for (i=0; i<NormalAxesCount; i++)		/* for all normal axes */
  {
    if (! NormalAxes[i]->order) continue;		/* an infinity axis */
    normalize (NormalAxes[i]->direction);
    for (k=0; k<DIMENSION; k++) norm[k] = NormalAxes[i]->direction[k];
    i_one = -1;
    n_one = n_zero = 0;
    for (k=0; k<DIMENSION; k++)
    {
      if (fabs(fabs(norm[k])-1) <= ToleranceFinal)
      {
	i_one = k;
	n_one ++;
      }
      if (fabs(norm[k]) <= ToleranceFinal) n_zero ++;
    }
    if (n_one!=1 || n_zero!=2)         continue;	/* not aligned axis */
    if (NormalAxes[i]->order != maxis) continue;   /* not principal axis */
    if (i_one != 2)		       continue;	/* not Z axis */

    /* Mark unique atoms preferring ones lying in the YZ plane */
    for (j=0; j<AtomsCount; j++)			/* for all atoms */
    {
      if (Atoms[j].type < 0)			continue;/*already non unique*/
      if (fabs(Atoms[j].x[0]) > ToleranceFinal) continue;    /* not YZ plane */
      rotate_atom (NormalAxes[i], Atoms+j, &atom);
      n_zero = 0;
      for (k=0; k<DIMENSION; k++)
	if (fabs(Atoms[j].x[k]-atom.x[k]) <= ToleranceFinal) n_zero ++;
      if (n_zero == 3) continue;      /* the original itself is close enough */
      for (m=1; m<NormalAxes[i]->order; m++) /*accumulate all rotation angles*/
      {
	for (n=0; n<AtomsCount; n++) /* search another close unmarked atom */
	{
	  if (n == j)		 continue;		/* original atom */
	  if (Atoms[n].type < 0) continue;   /* already marked as non unique */
	  n_zero = 0;
	  for (k=0; k<DIMENSION; k++)
	    if (fabs(Atoms[n].x[k]-atom.x[k]) <= ToleranceFinal) n_zero ++;
	  if (Atoms[n].type == atom.type) n_zero ++;	/* check atom number */
	  if (n_zero == 4)			/* this one is close enough */
	  {
	    Atoms[n].type = -Atoms[j].type;			/* mark it */
	    break;
	  }
	}
	rotate_atom (NormalAxes[i], &atom, &atom2);
	for (k=0; k<DIMENSION; k++) atom.x[k] = atom2.x[k];
      }
    }							/* for all atoms */

    /* Mark unique atoms preferring ones lying in the XZ plane */
    for (j=0; j<AtomsCount; j++)			/* for all atoms */
    {
      if (Atoms[j].type < 0)			continue;/*already non unique*/
      if (fabs(Atoms[j].x[1]) > ToleranceFinal) continue;    /* not XZ plane */
      rotate_atom (NormalAxes[i], Atoms+j, &atom);
      n_zero = 0;
      for (k=0; k<DIMENSION; k++)
	if (fabs(Atoms[j].x[k]-atom.x[k]) <= ToleranceFinal) n_zero ++;
      if (n_zero == 3) continue;      /* the original itself is close enough */
      for (m=1; m<NormalAxes[i]->order; m++) /*accumulate all rotation angles*/
      {
	for (n=0; n<AtomsCount; n++) /* search another close unmarked atom */
	{
	  if (n == j)		 continue;		/* original atom */
	  if (Atoms[n].type < 0) continue;   /* already marked as non unique */
	  n_zero = 0;
	  for (k=0; k<DIMENSION; k++)
	    if (fabs(Atoms[n].x[k]-atom.x[k]) <= ToleranceFinal) n_zero ++;
	  if (Atoms[n].type == atom.type) n_zero ++;	/* check atom number */
	  if (n_zero == 4)			/* this one is close enough */
	  {
	    Atoms[n].type = -Atoms[j].type;			/* mark it */
	    break;
	  }
	}
	rotate_atom (NormalAxes[i], &atom, &atom2);
	for (k=0; k<DIMENSION; k++) atom.x[k] = atom2.x[k];
      }
    }							/* for all atoms */
  }						/* for all normal axes */
}
  
/*************************************************************************
 * Mark atoms symmetric around either all or only abelian improper axes
 * to find out symmetry unique atoms.
 * For each atom in sequence other atoms are searched for which are
 * sufficiently close to all its images.
 * If they are close enough, they are marked.
 * There is some similarity with symmetrize_improper_axes()
 * but the criterium of closeness is different.
 *************************************************************************/
static void mark_symmetric_improper_axes (int check_abel)
{
  int	i, j, k, m, n, n_one, n_zero;
  ATOM1	atom, atom2;
  
  for (i=0; i<ImproperAxesCount; i++)		/* for all improper axes */
  {
    if (check_abel)			/* abelian improper axes only ? */
    {
      if (ImproperAxes[i]->order != 2) continue;     /* only order 2 allowed */
      if (fabs(ImproperAxes[i]->distance) > ToleranceFinal) continue; /* bad */
      n_one = n_zero = 0;
      for (k=0; k<DIMENSION; k++)
      {
	if (fabs(fabs(ImproperAxes[i]->direction[k])-1) <= ToleranceFinal)
	  n_one ++;
	if (fabs(ImproperAxes[i]->direction[k]) <= ToleranceFinal) n_zero ++;
      }
      if (n_one!=1 || n_zero!=2) continue;		/* not abelian */
    }

    for (j=0; j<AtomsCount; j++)			/* for all atoms */
    {
      if (Atoms[j].type < 0) continue;	/* already marked as non unique */
      rotate_reflect_atom (ImproperAxes[i], Atoms+j, &atom);
      n_zero = 0;
      for (k=0; k<DIMENSION; k++)
	if (fabs(Atoms[j].x[k]-atom.x[k]) <= ToleranceFinal) n_zero ++;
      if (n_zero == 3) continue;      /* the original itself is close enough */
      for (m=1; m<ImproperAxes[i]->order; m++)	/* accumulate all images */
      {
	for (n=j+1; n<AtomsCount; n++) /* search another close unmarked atom */
	{
	  if (Atoms[n].type < 0) continue;   /* already marked as non unique */
	  n_zero = 0;
	  for (k=0; k<DIMENSION; k++)
	    if (fabs(Atoms[n].x[k]-atom.x[k]) <= ToleranceFinal) n_zero ++;
	  if (Atoms[n].type == atom.type) n_zero ++;	/* check atom number */
	  if (n_zero == 4)			/* this one is close enough */
	  {
	    Atoms[n].type = -Atoms[j].type;			/* mark it */
	    break;
	  }
	}
	rotate_reflect_atom (ImproperAxes[i], &atom, &atom2);
	for (k=0; k<DIMENSION; k++) atom.x[k] = atom2.x[k];
      }
    }							/* for all atoms */
  }						/* for all improper axes */
}
  
/*************************************************************************
 * Mark atoms symmetric around either all or only abelian inversion centers
 * to find out symmetry unique atoms.
 * For each atom in sequence another atom is searched for which is
 * sufficiently close to its inversion. If it is close enough, it is marked.
 * There is some similarity with symmetrize_inversion_centers()
 * but the criterium of closeness is different.
 *************************************************************************/
static void mark_symmetric_inversion_centers (int check_abel)
{
  int	i, j, k, n, n_zero;
  ATOM1	atom;
  
  for (i=0; i<InversionCentersCount; i++)		/* for all centers */
  {
    if (check_abel && fabs(InversionCenters[i]->distance) > ToleranceFinal)
      continue;					/* abelian centers only ? */

    for (j=0; j<AtomsCount; j++)			/* for all atoms */
    {
      if (Atoms[j].type < 0) continue;	/* already marked as non unique */
      invert_atom (InversionCenters[i], Atoms+j, &atom);
      n_zero = 0;
      for (k=0; k<DIMENSION; k++)
	if (fabs(Atoms[j].x[k]-atom.x[k]) <= ToleranceFinal) n_zero ++;
      if (n_zero == 3) continue;      /* the original itself is close enough */
      for (n=j+1; n<AtomsCount; n++)   /* search another close unmarked atom */
      {
	if (Atoms[n].type < 0) continue;     /* already marked as non unique */
	n_zero = 0;
	for (k=0; k<DIMENSION; k++)
	  if (fabs(Atoms[n].x[k]-atom.x[k]) <= ToleranceFinal) n_zero ++;
	if (Atoms[n].type == atom.type) n_zero ++;	/* check atom number */
	if (n_zero == 4)			/* this one is close enough */
	{
	  Atoms[n].type = -Atoms[j].type;			/* mark it */
	  break;
	}
      }
    }							/* for all atoms */
  }							/* for all centers */
}

/*************************************************************************
 * Rank symmety planes according to the number of symmetry unique atoms.
 * For each atom in sequence another atom is searched for which is
 * sufficiently close to its mirror.
 * If it is close enough, it is marked and the rank counted down.
 * Very similar to mark_symmetric_planes().
 *************************************************************************/
static void rank_symmetric_planes (void)
{
  int	i, j, k, n, n_zero;
  ATOM1	atom;
  
  for (i=0; i<PlanesCount; i++)				/* for all planes */
  {
    Planes[i]->rank = AtomsCount;	/* initialize for all atoms unique */
    for (j=0; j<AtomsCount; j++) Atoms[j].type = abs(Atoms[j].type); /*unmark*/

    for (j=0; j<AtomsCount; j++)			/* for all atoms */
    {
      if (Atoms[j].type < 0) continue;	/* already marked as non unique */
      mirror_atom (Planes[i], Atoms+j, &atom);
      n_zero = 0;
      for (k=0; k<DIMENSION; k++)
	if (fabs(Atoms[j].x[k]-atom.x[k]) <= ToleranceFinal) n_zero ++;
      if (n_zero == 3) continue;      /* the original itself is close enough */
      for (n=j+1; n<AtomsCount; n++)   /* search another close unmarked atom */
      {
	if (Atoms[n].type < 0) continue;     /* already marked as non unique */
	n_zero = 0;
	for (k=0; k<DIMENSION; k++)
	  if (fabs(Atoms[n].x[k]-atom.x[k]) <= ToleranceFinal) n_zero ++;
	if (Atoms[n].type == atom.type) n_zero ++;	/* check atom number */
	if (n_zero == 4)			/* this one is close enough */
	{
	  Atoms[n].type = -Atoms[j].type;			/* mark it */
	  Planes[i]->rank --;			/* and count rank down */
	  break;
	}
      }
    }							/* for all atoms */
  }							/* for all planes */

  for (j=0; j<AtomsCount; j++) Atoms[j].type = abs(Atoms[j].type); /* unmark */
}

/*************************************************************************
 * Rank symmety axes according to the number of symmetry unique atoms.
 * For each atom in sequence other atoms are searched for which are
 * sufficiently close to all its rotation images.
 * If they are close enough, they are marked and the rank counted down.
 * Very similar to mark_symmetric_axes().
 *************************************************************************/
static void rank_symmetric_axes (void)
{
  int	i, j, k, m, n, n_zero;
  ATOM1	atom, atom2;
  
  for (i=0; i<NormalAxesCount; i++)		/* for all normal axes */
  {
    NormalAxes[i]->rank = AtomsCount;	/* initialize for all atoms unique */

    if (! NormalAxes[i]->order) continue;    /* nothing to mark for infinity */

    for (j=0; j<AtomsCount; j++) Atoms[j].type = abs(Atoms[j].type); /*unmark*/

    for (j=0; j<AtomsCount; j++)			/* for all atoms */
    {
      if (Atoms[j].type < 0) continue;	/* already marked as non unique */
      rotate_atom (NormalAxes[i], Atoms+j, &atom);
      n_zero = 0;
      for (k=0; k<DIMENSION; k++)
	if (fabs(Atoms[j].x[k]-atom.x[k]) <= ToleranceFinal) n_zero ++;
      if (n_zero == 3) continue;      /* the original itself is close enough */
      for (m=1; m<NormalAxes[i]->order; m++) /*accumulate all rotation angles*/
      {
	for (n=j+1; n<AtomsCount; n++) /* search another close unmarked atom */
	{
	  if (Atoms[n].type < 0) continue;   /* already marked as non unique */
	  n_zero = 0;
	  for (k=0; k<DIMENSION; k++)
	    if (fabs(Atoms[n].x[k]-atom.x[k]) <= ToleranceFinal) n_zero ++;
	  if (Atoms[n].type == atom.type) n_zero ++;	/* check atom number */
	  if (n_zero == 4)			/* this one is close enough */
	  {
	    Atoms[n].type = -Atoms[j].type;			/* mark it */
	    NormalAxes[i]->rank --;		/* and count rank down */
	    break;
	  }
	}
	rotate_atom (NormalAxes[i], &atom, &atom2);
	for (k=0; k<DIMENSION; k++) atom.x[k] = atom2.x[k];
      }
    }							/* for all atoms */
  }						/* for all normal axes */

  for (j=0; j<AtomsCount; j++) Atoms[j].type = abs(Atoms[j].type); /* unmark */
}
  
/*************************************************************************
 * Rank improper axes according to the number of symmetry unique atoms.
 * For each atom in sequence other atoms are searched for which are
 * sufficiently close to all its images.
 * If they are close enough, they are marked and the rank counted down.
 * Very similar to mark_symmetric_improper_axes().
 *************************************************************************/
static void rank_symmetric_improper_axes (void)
{
  int	i, j, k, m, n, n_zero;
  ATOM1	atom, atom2;
  
  for (i=0; i<ImproperAxesCount; i++)		/* for all improper axes */
  {
    ImproperAxes[i]->rank = AtomsCount;	/* initialize for all atoms unique */
    for (j=0; j<AtomsCount; j++) Atoms[j].type = abs(Atoms[j].type); /*unmark*/

    for (j=0; j<AtomsCount; j++)			/* for all atoms */
    {
      if (Atoms[j].type < 0) continue;	/* already marked as non unique */
      rotate_reflect_atom (ImproperAxes[i], Atoms+j, &atom);
      n_zero = 0;
      for (k=0; k<DIMENSION; k++)
	if (fabs(Atoms[j].x[k]-atom.x[k]) <= ToleranceFinal) n_zero ++;
      if (n_zero == 3) continue;      /* the original itself is close enough */
      for (m=1; m<ImproperAxes[i]->order; m++)	/* accumulate all images */
      {
	for (n=j+1; n<AtomsCount; n++) /* search another close unmarked atom */
	{
	  if (Atoms[n].type < 0) continue;   /* already marked as non unique */
	  n_zero = 0;
	  for (k=0; k<DIMENSION; k++)
	    if (fabs(Atoms[n].x[k]-atom.x[k]) <= ToleranceFinal) n_zero ++;
	  if (Atoms[n].type == atom.type) n_zero ++;	/* check atom number */
	  if (n_zero == 4)			/* this one is close enough */
	  {
	    Atoms[n].type = -Atoms[j].type;			/* mark it */
	    ImproperAxes[i]->rank --;		/* and count rank down */
	    break;
	  }
	}
	rotate_reflect_atom (ImproperAxes[i], &atom, &atom2);
	for (k=0; k<DIMENSION; k++) atom.x[k] = atom2.x[k];
      }
    }							/* for all atoms */
  }						/* for all improper axes */

  for (j=0; j<AtomsCount; j++) Atoms[j].type = abs(Atoms[j].type); /* unmark */
}
  
/*************************************************************************
 * Rank inversion centers according to the number of symmetry unique atoms.
 * For each atom in sequence another atom is searched for which is
 * sufficiently close to its inversion.
 * If it is close enough, it is marked and the rank counted down.
 * Very similar to mark_symmetric_inversion_centers().
 *************************************************************************/
static void rank_symmetric_inversion_centers (void)
{
  int	i, j, k, n, n_zero;
  ATOM1	atom;
  
  for (i=0; i<InversionCentersCount; i++)		/* for all centers */
  {
    InversionCenters[i]->rank = AtomsCount; /*initialize for all atoms unique*/
    for (j=0; j<AtomsCount; j++) Atoms[j].type = abs(Atoms[j].type); /*unmark*/

    for (j=0; j<AtomsCount; j++)			/* for all atoms */
    {
      if (Atoms[j].type < 0) continue;	/* already marked as non unique */
      invert_atom (InversionCenters[i], Atoms+j, &atom);
      n_zero = 0;
      for (k=0; k<DIMENSION; k++)
	if (fabs(Atoms[j].x[k]-atom.x[k]) <= ToleranceFinal) n_zero ++;
      if (n_zero == 3) continue;      /* the original itself is close enough */
      for (n=j+1; n<AtomsCount; n++)   /* search another close unmarked atom */
      {
	if (Atoms[n].type < 0) continue;     /* already marked as non unique */
	n_zero = 0;
	for (k=0; k<DIMENSION; k++)
	  if (fabs(Atoms[n].x[k]-atom.x[k]) <= ToleranceFinal) n_zero ++;
	if (Atoms[n].type == atom.type) n_zero ++;	/* check atom number */
	if (n_zero == 4)			/* this one is close enough */
	{
	  Atoms[n].type = -Atoms[j].type;			/* mark it */
	  InversionCenters[i]->rank --;		/* and count rank down */
	  break;
	}
      }
    }							/* for all atoms */
  }							/* for all centers */

  for (j=0; j<AtomsCount; j++) Atoms[j].type = abs(Atoms[j].type); /* unmark */
}

/*************************************************************************
 * Canonize atomic numbers and coordinates against symmetry unique flags
 * and some round off errors near zero
 *************************************************************************/
static void canonize_atoms (void)
{
  int		i, k;
  double	prec;

  prec = ToleranceFinal;
  if (prec > .001) prec = .001;		/* minimum acceptable precision */
  for (i=0; i<AtomsCount; i++)
  {
    Atoms[i].type = abs(Atoms[i].type);	/* unmark symmetry flag */
    for (k=0; k<DIMENSION; k++)	/* make almost zero coordinate zero exactly */
      if (fabs(Atoms[i].x[k]) < prec) Atoms[i].x[k] = 0;
  }
}

/*************************************************************************
 * Print number of and coordinates of all symmetry unique (not marked) atoms
 *************************************************************************/
static void report_all_atoms (void)
{
  int i, n, width, prec;

  prec = (int)ceil(-log10(ToleranceFinal));   /* derive reasonable precision */
  if (prec < 3) prec = 3;		/* minimum acceptable precision */
  width = prec+5;
  n = 0;
  for (i=0; i<AtomsCount; i++) if (Atoms[i].type >= 0) n++;   /* no of atoms */
  printf ("%d\n", n);
  for (i=0; i<AtomsCount; i++)
    if (Atoms[i].type >= 0)				/* symmetry unique ? */
      printf ("%4d %*.*f %*.*f %*.*f\n", Atoms[i].type,
	      width, prec, Atoms[i].x[0],
	      width, prec, Atoms[i].x[1],
	      width, prec, Atoms[i].x[2]);
}

/*************************************************************************
 * Print abelian symmetry planes. Derived from report_planes().
 *************************************************************************/
static void report_abel_planes (void)
{
  int i, k, n_zero, n_one, nplanes;

  nplanes = 0;
  for (i=0; i<PlanesCount; i++)				/* count planes */
  {
    if (fabs(Planes[i]->distance) > ToleranceFinal) continue; /* not abelian */
    n_one = n_zero = 0;
    for (k=0; k<DIMENSION; k++)
    {
      if (fabs(fabs(Planes[i]->normal[k])-1) <= ToleranceFinal) n_one  ++;
      if (fabs(Planes[i]->normal[k])         <= ToleranceFinal) n_zero ++;
    }
    if (n_one==1 && n_zero==2) nplanes ++;
  }

  if (nplanes == 0)
    printf ("There are no abelian planes of symmetry in the molecule\n");
  else
  {
    if (nplanes == 1)
      printf ("There is an abelian plane of symmetry in the molecule\n");
    else printf ("There are %d abelian planes of symmetry in the molecule\n",
		 nplanes);
    printf ("    Rank Residual          Direction of the normal           Distance\n");
    for (i=0; i<PlanesCount; i++)			/* print planes */
    {
      n_one = n_zero = 0;
      for (k=0; k<DIMENSION; k++)
      {
	if (fabs(fabs(Planes[i]->normal[k])-1) <= ToleranceFinal) n_one  ++;
	if (fabs(Planes[i]->normal[k])         <= ToleranceFinal) n_zero ++;
      }
      if (n_one!=1 || n_zero!=2) continue;		/* not abelian */
      printf ("%3d%4d %8.4e ", i+1, Planes[i]->rank, Planes[i]->maxdev);
      printf ("(%11.8f,%11.8f,%11.8f) ",
	      Planes[i]->normal[0],
	      Planes[i]->normal[1],
	      Planes[i]->normal[2]);
      printf ("%14.8f\n", Planes[i]->distance);
    }
  }
}

/*************************************************************************
 * Print abelian normal axes. Derived from report_axes().
 *************************************************************************/
static void report_abel_axes (void)
{
  int i, k, n_zero, n_one, naxes;

  naxes = 0;
  for (i=0; i<NormalAxesCount; i++)			/* count normal axes */
  {					/* only order 2 or inf allowed */
    if (NormalAxes[i]->order && NormalAxes[i]->order != 2) continue;
    if (fabs(NormalAxes[i]->distance) > ToleranceFinal)    continue;  /* bad */
    n_one = n_zero = 0;
    for (k=0; k<DIMENSION; k++)
    {
      if (fabs(fabs(NormalAxes[i]->direction[k])-1) <= ToleranceFinal)
	n_one ++;
      if (fabs(NormalAxes[i]->direction[k]) <= ToleranceFinal) n_zero ++;
    }
    if (n_one==1 && n_zero==2) naxes ++;
  }

  if (naxes == 0)
    printf ("There are no abelian normal axes in the molecule\n");
  else
  {
    if (naxes == 1)
      printf ("There is an abelian normal axis in the molecule\n");
    else
      printf ("There are %d abelian normal axes in the molecule\n", naxes);
    printf ("    Rank Residual  Order         Direction of the axis                         Supporting point\n");
    for (i=0; i<NormalAxesCount; i++)			/* print normal axes */
    {					/* only order 2 or inf allowed */
      if (NormalAxes[i]->order && NormalAxes[i]->order != 2) continue;
      n_one = n_zero = 0;
      for (k=0; k<DIMENSION; k++)
      {
	if (fabs(fabs(NormalAxes[i]->direction[k])-1) <= ToleranceFinal)
	  n_one ++;
	if (fabs(NormalAxes[i]->direction[k]) <= ToleranceFinal) n_zero ++;
      }
      if (n_one!=1 || n_zero!=2) continue;		/* not abelian */
      printf ("%3d%4d %8.4e ",
	      i+1, NormalAxes[i]->rank, NormalAxes[i]->maxdev);
      if (NormalAxes[i]->order)	printf ("%3d ", NormalAxes[i]->order);
      else printf ("Inf ");
      printf ("(%11.8f,%11.8f,%11.8f) ",
	      NormalAxes[i]->direction[0],
	      NormalAxes[i]->direction[1],
	      NormalAxes[i]->direction[2]);
      printf ("(%14.8f,%14.8f,%14.8f)\n",
	      NormalAxes[0]->distance*NormalAxes[0]->normal[0],
	      NormalAxes[0]->distance*NormalAxes[0]->normal[1],
	      NormalAxes[0]->distance*NormalAxes[0]->normal[2]);
    }
  }
}

/*************************************************************************
 * Print abelian improper axes. Derived from report_improper_axes().
 *************************************************************************/
static void report_abel_improper_axes (void)
{
  int i, k, n_zero, n_one, naxes;

  naxes = 0;
  for (i=0; i<ImproperAxesCount; i++)		/* count improper axes */
  {						/* only order 2 allowed */
    if (ImproperAxes[i]->order != 2)                      continue;
    if (fabs(ImproperAxes[i]->distance) > ToleranceFinal) continue;   /* bad */
    n_one = n_zero = 0;
    for (k=0; k<DIMENSION; k++)
    {
      if (fabs(fabs(ImproperAxes[i]->direction[k])-1) <= ToleranceFinal)
	n_one ++;
      if (fabs(ImproperAxes[i]->direction[k]) <= ToleranceFinal) n_zero ++;
    }
    if (n_one==1 && n_zero==2) naxes ++;
  }

  if (naxes == 0)
    printf ("There are no abelian improper axes in the molecule\n");
  else
  {
    if (naxes == 1)
      printf ("There is an abelian improper axis in the molecule\n");
    else
      printf ("There are %d abelian improper axes in the molecule\n", naxes);
    printf ("    Rank Residual  Order         Direction of the axis                         Supporting point\n");
    for (i=0; i<ImproperAxesCount; i++)		/* print improper axes */
    {
      if (ImproperAxes[i]->order != 2) continue;     /* only order 2 allowed */
      n_one = n_zero = 0;
      for (k=0; k<DIMENSION; k++)
      {
	if (fabs(fabs(ImproperAxes[i]->direction[k])-1) <= ToleranceFinal)
	  n_one ++;
	if (fabs(ImproperAxes[i]->direction[k]) <= ToleranceFinal) n_zero ++;
      }
      if (n_one!=1 || n_zero!=2) continue;		/* not abelian */
      printf ("%3d%4d %8.4e ",
	      i+1, ImproperAxes[i]->rank, ImproperAxes[i]->maxdev);
      printf ("%3d ", ImproperAxes[i]->order);
      printf ("(%11.8f,%11.8f,%11.8f) ",
	      ImproperAxes[i]->direction[0],
	      ImproperAxes[i]->direction[1],
	      ImproperAxes[i]->direction[2]);
      printf ("(%14.8f,%14.8f,%14.8f)\n",
	      ImproperAxes[0]->distance*ImproperAxes[0]->normal[0],
	      ImproperAxes[0]->distance*ImproperAxes[0]->normal[1],
	      ImproperAxes[0]->distance*ImproperAxes[0]->normal[2]);
    }
  }
}

/*************************************************************************
 * Print abelian inversion centers. Derived from report_inversion_centers().
 *************************************************************************/
static void report_abel_inversion_centers (void)
{
  if (InversionCentersCount == 0 ||
      fabs(InversionCenters[0]->distance) > ToleranceFinal)
    printf ("There is no abelian inversion center in the molecule\n");
  else
  {
    printf ("There is an abelian inversion center in the molecule\n");
    printf (" Rank Residual                      Position\n");
    printf ("%4d %8.4e ",
	    InversionCenters[0]->rank, InversionCenters[0]->maxdev);
    printf ("(%14.8f,%14.8f,%14.8f)\n",
	    InversionCenters[0]->distance*InversionCenters[0]->normal[0],
	    InversionCenters[0]->distance*InversionCenters[0]->normal[1],
	    InversionCenters[0]->distance*InversionCenters[0]->normal[2]);
  }
}

/*************************************************************************
 * Print abelian symmetry elements, symmetry generators (IGLO style)
 * and determine abelian symmetry subgroup
 * Derived from report_symmetry_elements_brief()
 *************************************************************************/
static void report_abel_symmetry_elements_brief (void)
{
  int	i, j, k, n, n_zero, n_one, i_one,
	nplanes, naxes, nimproper, ninf, ncenters;
  int	splane[DIMENSION], saxe[DIMENSION], simproper[DIMENSION],
	sinf[DIMENSION];
  char	buf[100];
  char	*symmetry_code;
  char	*genplane[DIMENSION] = {  "x",  "y",  "z" };
  char	*genaxe  [DIMENSION] = { "yz", "xz", "xy" };

  nplanes = naxes = nimproper = ncenters = ninf = 0;
  for (k=0; k<DIMENSION; k++) splane[k] = saxe[k] = simproper[k] = sinf[k] = 0;

  /* Count symmetry elements... */
  for (i=0; i<PlanesCount; i++)			/* count abelian planes */
  {
    if (fabs(Planes[i]->distance) > ToleranceFinal) continue; /* not abelian */
    i_one = -1;
    n_one = n_zero = 0;
    for (k=0; k<DIMENSION; k++)
    {
      if (fabs(fabs(Planes[i]->normal[k])-1) <= ToleranceFinal)
      {
	i_one = k;
	n_one ++;
      }
      if (fabs(Planes[i]->normal[k]) <= ToleranceFinal) n_zero ++;
    }
    if (n_one!=1 || n_zero!=2) continue;		/* not abelian */
    nplanes ++;
    splane[i_one] = 1;				/* mark symmetry generator */
  }

  for (i=0; i<NormalAxesCount; i++)		/* count abelian normal axes */
  {					/* only order 2 or inf allowed */
    if (NormalAxes[i]->order && NormalAxes[i]->order != 2) continue;
    if (fabs(NormalAxes[i]->distance) > ToleranceFinal)    continue;  /* bad */
    i_one = -1;
    n_one = n_zero = 0;
    for (k=0; k<DIMENSION; k++)
    {
      if (fabs(fabs(NormalAxes[i]->direction[k])-1) <= ToleranceFinal)
      {
	i_one = k;
	n_one ++;
      }
      if (fabs(NormalAxes[i]->direction[k]) <= ToleranceFinal) n_zero ++;
    }
    if (n_one!=1 || n_zero!=2) continue;		/* not abelian */
    if (NormalAxes[i]->order)
    {
      naxes ++;
      saxe[i_one] = 1;			/* mark normal symmetry generator */
    }
    else
    {
      ninf ++;
      sinf[i_one] = 1;			/* mark infinity symmetry generator */
    }
  }

  for (i=0; i<ImproperAxesCount; i++)	/* count abelian improper axes */
  {						/* only order 2 allowed */
    if (ImproperAxes[i]->order != 2)                      continue;
    if (fabs(ImproperAxes[i]->distance) > ToleranceFinal) continue;   /* bad */
    i_one = -1;
    n_one = n_zero = 0;
    for (k=0; k<DIMENSION; k++)
    {
      if (fabs(fabs(ImproperAxes[i]->direction[k])-1) <= ToleranceFinal)
      {
	i_one = k;
	n_one ++;
      }
      if (fabs(ImproperAxes[i]->direction[k]) <= ToleranceFinal) n_zero ++;
    }
    if (n_one!=1 || n_zero!=2) continue;		/* not abelian */
    nimproper ++;
    simproper[i_one] = 1;			/* mark symmetry generator */
  }

  for (i=0; i<InversionCentersCount; i++) /* count abelian inversion centers */
    if (fabs(InversionCenters[i]->distance) <= ToleranceFinal) ncenters ++;

  /* Derive corresponding symmetry elements and generators... */
  for (k=0; k<DIMENSION; k++)
    if (sinf[k])     /* there are always two (Cinfv) or three (Dinfh) planes */
    {
      for (j=0; j<DIMENSION; j++) if (j != k || ncenters) splane[j] = 1;
      sinf[k] = 0;			/* clear unneeded infinity axes */
    }
  for (k=0; k<DIMENSION; k++)
    if (splane[k])				/* clear extra unneeded axes */
    {
      for (j=0; j<DIMENSION; j++) if (j != k) saxe[j] = 0;
      ncenters = 0;			/* clear unneeded inversion centers */
    }
  for (k=0; k<DIMENSION; k++) if (saxe[k]) ncenters = 0;    /* clear centers */

  /* Compose symmetry code... */
  symmetry_code =
    calloc (1, 10*(nplanes+naxes+nimproper+ninf+InversionCentersCount+2));
  if (symmetry_code == NULL)
  {
    fprintf (stderr, "Unable to allocate memory for abelian symmetry ID code in report_abel_symmetry_elements_brief()\n");
    exit (EXIT_FAILURE);
  }
  if (nplanes+naxes+nimproper+ninf+InversionCentersCount == 0)
    printf ("Molecule has no abelian symmetry elements\n");
  else
  {
    printf ("Molecule has the following abelian symmetry elements: ");
    if (InversionCentersCount > 0) strcat (symmetry_code, "(i) ");
    if (ninf == 1)                 strcat (symmetry_code, "(Cinf) ");
    if (ninf > 1)
    {
      sprintf (buf, "%d*(Cinf) ", ninf);
      strcat  (symmetry_code, buf);
    }
    if (naxes == 1)
    {
      sprintf (buf, "(C2) ");
      strcat  (symmetry_code, buf);
    }
    if (naxes > 1)
    {
      sprintf (buf, "%d*(C2) ", naxes);
      strcat  (symmetry_code, buf);
    }
    if (nimproper == 1)
    {
      sprintf (buf, "(S2) ");
      strcat  (symmetry_code, buf);
    }
    if (nimproper > 1)
    {
      sprintf (buf, "%d*(S2) ", nimproper);
      strcat  (symmetry_code, buf);
    }
    if (nplanes == 1) strcat (symmetry_code, "(sigma) ");
    if (nplanes > 1)
    {
      sprintf (buf, "%d*(sigma) ", nplanes);
      strcat  (symmetry_code, buf);
    }
    printf ("%s\n", symmetry_code);

    /* Print symmetry generators */
    printf ("Molecule has the following abelian symmetry generators: ");
    n = 0;
    for (k=0; k<DIMENSION; k++)
    {
      if (splane[k])
      {
	if (n) printf (",");
	printf ("%s", genplane[k]);
	n ++;
      }
      if (saxe[DIMENSION-k-1])
      {
	if (n) printf (",");
	printf ("%s", genaxe[DIMENSION-k-1]);
	n ++;
      }
    }
    if (ncenters)
    {
      if (n) printf (",");
      printf ("xyz");
      n ++;
    }
    printf ("\n");
  }

  SymmetryCode = symmetry_code;					/* Whew...! */
}

int point_group_calc(int num_atom, double *charge, double *coord, char *PG, double fin_tol){

    int i;
	
	AtomsCount = num_atom;
	
	ToleranceFinal = fin_tol;
	
    Atoms = calloc( num_atom, sizeof( ATOM1 ) ) ;
    if( Atoms == NULL ){
        fprintf( stderr, "Out of memory for atoms coordinates\n" ) ;
        return -1 ;
    }
    for( i = 0 ; i < num_atom ; i++ ){
	    Atoms[i].type = (int)charge[i];
	    Atoms[i].x[0] = coord[3*i];
	    Atoms[i].x[1] = coord[3*i+1];
	    Atoms[i].x[2] = coord[3*i+2];
    }
	
	find_symmetry_elements() ;
    sort_symmetry_elements() ;
    summarize_symmetry_elements() ;
    if( BadOptimization )
        printf( "Refinement of some symmetry elements was terminated before convergence was reached.\n"
                "Some symmetry elements may remain unidentified.\n" ) ;
    
    rank_symmetric_planes ();		/* particular order not important */
    rank_symmetric_axes ();
    rank_symmetric_inversion_centers ();
    rank_symmetric_improper_axes ();
    
    if( verbose >= 0 )
        report_symmetry_elements_verbose() ;
    report_symmetry_elements_brief() ;
    identify_point_group1( PG ) ;

	return 0;
}

/*************************************************************************
 * End SGE code
 *************************************************************************/

//int
//main( int argc, char **argv )
//{
//        char          *program = *argv ;
//        FILE          *in ;
//
//for( argc--, argv++ ; argc > 0 ; argc -= 2, argv += 2 ){
//    if( **argv != '-' )
//        break ;
//    if( strcmp( *argv, "-help"         ) == 0 ||
//        strcmp( *argv, "-h"            ) == 0 ||
//        strcmp( *argv, "-?"            ) == 0 ){
//        argc++ ; argv-- ;
//        printf( "%s [option value ...] [filename]\n" 
//                "Valid options are:\n"
//                "  -verbose      (%3d) Determines verbosity level\n"
//                "                      All values above 0 are intended for debugging purposes\n"
//                "  -maxaxisorder (%3d) Maximum order of rotation axis to look for\n"
//                "  -maxoptcycles (%3d) Maximum allowed number of cycles in symmetry element optimization\n"
//                "  --                  Terminates option processing\n"
//                "Defaults should be Ok for these:\n"
//                "  -same         (%8g) Atoms are colliding if distance falls below this value\n"
//                "  -primary      (%8g) Initial loose criterion for atom equivalence\n"
//                "  -final        (%8g) Final criterion for atom equivalence\n"
//                "  -maxoptstep   (%8g) Largest step allowed in symmetry element optimization\n"
//                "  -minoptstep   (%8g) Termination criterion in symmetry element optimization\n"
//                "  -gradstep     (%8g) Finite step used in numeric gradient evaluation\n" 
//                "  -minchange    (%8g) Minimum allowed change in target function\n" 
//                "  -minchgcycles (%8d) Number of minchange cycles before optimization stops\n"
//                "  -maxuniq                 Reorient molecule preferring more unique atoms\n"
//                "  -minuniq      (default)  Reorient molecule preferring less unique atoms\n"
//                "  -symmplane      (all)    Symmetrize atoms according to this plane number only\n"
//                "  -symmaxis       (all)    Symmetrize atoms according to this axis number only\n"
//                "  -symmimproper   (all)    Symmetrize atoms according to this improper axis number only\n"
//                "  -symmcenter     (all)    Symmetrize atoms according to this inversion center number only\n"
//                "  -na                      Prefer non-abelian symmetry elements for reorientation\n"
//                "  -minaxis                 Prefer lower order axes for reorientation\n"
//                "                           (use together with -na for tetrahedral groups)\n",
//            program, verbose, MaxAxisOrder, MaxOptCycles, ToleranceSame, TolerancePrimary,
//            ToleranceFinal, MaxOptStep, MinOptStep, GradientStep, OptChangeThreshold, OptChangeHits ) ;
//        printf( "\n"
//                "Input is expected in the following format:\n"
//                "number_of_atoms\n"
//                "AtomicNumber X Y Z\n"
//                "...\n" ) ;
//        printf( "\n"
//                "Note that only primitive rotations will be reported\n" ) ;
//        printf( "This is version $Revision: 1.15 $ ($Date: 2000/01/25 16:47:17 $)\n" ) ;
//        exit( EXIT_SUCCESS ) ;
//        }
//    else
//    if( strcmp( *argv, "--"            ) == 0 ){
//        argc-- ; argv++ ; break ;
//        }
//    if( argc < 2 ){
//        fprintf( stderr, "Missing argument for \"%s\"\n", *argv ) ;
//        exit( EXIT_FAILURE ) ;
//        }
//    if( strcmp( *argv, "-minchgcycles" ) == 0 ){
//        if( sscanf( argv[1], "%d", &OptChangeHits ) != 1 ){
//            fprintf( stderr, "Invalid parameter for -minchgcycles: \"%s\"\n", argv[1] ) ;
//            exit( EXIT_FAILURE ) ;
//            }
//        }
//    else
//    if( strcmp( *argv, "-minchange"    ) == 0 ){
//        if( sscanf( argv[1], "%lg", &OptChangeThreshold ) != 1 ){
//            fprintf( stderr, "Invalid parameter for -minchange: \"%s\"\n", argv[1] ) ;
//            exit( EXIT_FAILURE ) ;
//            }
//        }
//    else
//    if( strcmp( *argv, "-same"         ) == 0 ){
//        if( sscanf( argv[1], "%lg", &ToleranceSame ) != 1 ){
//            fprintf( stderr, "Invalid parameter for -same: \"%s\"\n", argv[1] ) ;
//            exit( EXIT_FAILURE ) ;
//            }
//        }
//    else
//    if( strcmp( *argv, "-primary"      ) == 0 ){
//        if( sscanf( argv[1], "%lg", &TolerancePrimary ) != 1 ){
//            fprintf( stderr, "Invalid parameter for -primary: \"%s\"\n", argv[1] ) ;
//            exit( EXIT_FAILURE ) ;
//            }
//        }
//    else
//    if( strcmp( *argv, "-final"        ) == 0 ){
//        if( sscanf( argv[1], "%lg", &ToleranceFinal ) != 1 ){
//            fprintf( stderr, "Invalid parameter for -final: \"%s\"\n", argv[1] ) ;
//            exit( EXIT_FAILURE ) ;
//            }
//        }
//    else
//    if( strcmp( *argv, "-maxoptstep"   ) == 0 ){
//        if( sscanf( argv[1], "%lg", &MaxOptStep ) != 1 ){
//            fprintf( stderr, "Invalid parameter for -maxoptstep: \"%s\"\n", argv[1] ) ;
//            exit( EXIT_FAILURE ) ;
//            }
//        }
//    else
//    if( strcmp( *argv, "-minoptstep"   ) == 0 ){
//        if( sscanf( argv[1], "%lg", &MinOptStep ) != 1 ){
//            fprintf( stderr, "Invalid parameter for -minoptstep: \"%s\"\n", argv[1] ) ;
//            exit( EXIT_FAILURE ) ;
//            }
//        }
//    else
//    if( strcmp( *argv, "-gradstep"     ) == 0 ){
//        if( sscanf( argv[1], "%lg", &GradientStep ) != 1 ){
//            fprintf( stderr, "Invalid parameter for -gradstep: \"%s\"\n", argv[1] ) ;
//            exit( EXIT_FAILURE ) ;
//            }
//        }
//    else
//    if( strcmp( *argv, "-verbose"      ) == 0 ){
//        if( sscanf( argv[1], "%d", &verbose ) != 1 ){
//            fprintf( stderr, "Invalid parameter for -verbose: \"%s\"\n", argv[1] ) ;
//            exit( EXIT_FAILURE ) ;
//            }
//        }
//    else
//    if( strcmp( *argv, "-maxoptcycles" ) == 0 ){
//        if( sscanf( argv[1], "%d", &MaxOptCycles ) != 1 ){
//            fprintf( stderr, "Invalid parameter for -maxoptcycles: \"%s\"\n", argv[1] ) ;
//            exit( EXIT_FAILURE ) ;
//            }
//        }
//    else
//    if( strcmp( *argv, "-maxaxisorder" ) == 0 ){
//        if( sscanf( argv[1], "%d", &MaxAxisOrder ) != 1 ){
//            fprintf( stderr, "Invalid parameter for -maxaxisorder: \"%s\"\n", argv[1] ) ;
//            exit( EXIT_FAILURE ) ;
//            }
//        }
//    else
//    if( strcmp( *argv, "-maxuniq" ) == 0 ){
//	OrdUniq = -1;
//        argc++ ; argv-- ;
//        }
//    else
//    if( strcmp( *argv, "-minuniq" ) == 0 ){
//	OrdUniq = 1;
//        argc++ ; argv-- ;
//        }
//    else
//    if( strcmp( *argv, "-minaxis" ) == 0 ){
//	OrdAxis = -1;
//        argc++ ; argv-- ;
//        }
//    else
//    if( strcmp( *argv, "-symmplane" ) == 0 ){
//        if( sscanf( argv[1], "%d", &SymmPlane ) != 1 ){
//            fprintf( stderr, "Invalid parameter for -symmplane: \"%s\"\n", argv[1] ) ;
//            exit( EXIT_FAILURE ) ;
//            }
//	else SymmAxis = SymmImproper = SymmCenter = 0;
//        }
//    else
//    if( strcmp( *argv, "-symmaxis" ) == 0 ){
//        if( sscanf( argv[1], "%d", &SymmAxis ) != 1 ){
//            fprintf( stderr, "Invalid parameter for -symmaxis: \"%s\"\n", argv[1] ) ;
//            exit( EXIT_FAILURE ) ;
//            }
//	else SymmPlane = SymmImproper = SymmCenter = 0;
//        }
//    else
//    if( strcmp( *argv, "-symmimproper" ) == 0 ){
//        if( sscanf( argv[1], "%d", &SymmImproper ) != 1 ){
//            fprintf( stderr, "Invalid parameter for -symmimproper: \"%s\"\n", argv[1] ) ;
//            exit( EXIT_FAILURE ) ;
//            }
//	else SymmPlane = SymmAxis = SymmCenter = 0;
//        }
//    else
//    if( strcmp( *argv, "-symmcenter" ) == 0 ){
//        if( sscanf( argv[1], "%d", &SymmCenter ) != 1 ){
//            fprintf( stderr, "Invalid parameter for -symmcenter: \"%s\"\n", argv[1] ) ;
//            exit( EXIT_FAILURE ) ;
//            }
//	else SymmPlane = SymmAxis = SymmImproper = 0;
//        }
//    else
//    if( strcmp( *argv, "-na" ) == 0 ){
//	NonAbel = 1;
//        argc++ ; argv-- ;
//        }
//    else {
//        fprintf( stderr, "Unrecognized option \"%s\"\n", *argv ) ;
//        exit( EXIT_FAILURE ) ;
//        }
//    }
//if( argc > 0 ){
//    if( ( in = fopen( *argv, "rt" ) ) == NULL ){
//        perror( *argv ) ;
//        exit( EXIT_FAILURE ) ;
//        }
//    }
//else {
//    in = stdin ;
//    }
//if( read_coordinates( in ) < 0 ){
//    fprintf( stderr, "Error reading in atomic coordinates\n" ) ;
//    exit( EXIT_FAILURE ) ;
//    }
//fclose( in ) ;
//find_symmetry_elements() ;
//sort_symmetry_elements() ;
//summarize_symmetry_elements() ;
//if( BadOptimization )
//    printf( "Refinement of some symmetry elements was terminated before convergence was reached.\n"
//            "Some symmetry elements may remain unidentified.\n" ) ;
//
//rank_symmetric_planes ();		/* particular order not important */
//rank_symmetric_axes ();
//rank_symmetric_inversion_centers ();
//rank_symmetric_improper_axes ();
//
//if( verbose >= 0 )
//    report_symmetry_elements_verbose() ;
//report_symmetry_elements_brief() ;
//identify_point_group() ;
//
///*************************************************************************
// * Begin SGE code
// *************************************************************************/
//
//printf ("-----------------------------\n");
//printf ("Original coordinates of atoms\n");
//printf ("-----------------------------\n");
//report_all_atoms ();
//
//if (NonAbel)				/* abelian properties uninteresting */
//{
//  printf ("====================================================\n"); 
//  printf ("Reorienting the molecule along its symmetry elements\n"); 
//  printf ("====================================================\n"); 
//  reorient_planes ();	/* best in this particular order as planes preferred */
//  reorient_axes (0);
//  reorient_inversion_centers ();
//  reorient_improper_axes ();
//  canonize_axes (0);
//  if (verbose >= 0) report_symmetry_elements_verbose ();
//  printf ("-------------------------------\n");
//  printf ("Reoriented coordinates of atoms\n");
//  printf ("-------------------------------\n");
//  report_all_atoms ();
//  symmetrize_axes (0);	/* not abelian but all symmetry elements used */
//  symmetrize_improper_axes (0);
//  symmetrize_inversion_centers (0);
//  symmetrize_planes (0);  /* best in this reversed order as planes preferred */
//  canonize_atoms (); 
//  printf ("------------------------------------\n");
//  printf ("Symmetrized coordinates of all atoms\n");
//  printf ("------------------------------------\n");
//  report_all_atoms ();
//  mark_principal_axes ();			/* this should come first */
//  mark_symmetric_planes (0);		/* particular order not important */
//  mark_symmetric_axes (0);
//  mark_symmetric_inversion_centers (0);
//  mark_symmetric_improper_axes (0);
//  printf ("------------------------------------------\n");
//  printf ("Coordinates of total symmetry unique atoms\n");
//  printf ("------------------------------------------\n");
//  report_all_atoms ();
//}
//else					/* abelian properties requested */
//{
//  printf ("====================================================\n"); 
//  printf ("Reorienting the molecule along its symmetry elements\n"); 
//  printf ("====================================================\n"); 
//  reorient_planes ();	/* best in this particular order as planes preferred */
//  reorient_axes (1);
//  reorient_inversion_centers ();
//  reorient_improper_axes ();
//  canonize_axes (1);
//  if (verbose >= 0) report_symmetry_elements_verbose ();
//  printf ("-------------------------------\n");
//  printf ("Reoriented coordinates of atoms\n");
//  printf ("-------------------------------\n");
//  report_all_atoms ();
//  symmetrize_axes (0);	/* not abelian but all symmetry elements used */
//  symmetrize_improper_axes (0);
//  symmetrize_inversion_centers (0);
//  symmetrize_planes (0);  /* best in this reversed order as planes preferred */
//  canonize_atoms (); 
//  printf ("------------------------------------\n");
//  printf ("Symmetrized coordinates of all atoms\n");
//  printf ("------------------------------------\n");
//  report_all_atoms ();
//
//  printf ("==================================\n"); 
//  printf ("Abelian symmetry subgroup analysis\n"); 
//  printf ("==================================\n"); 
//  if (verbose >= 0)
//  {
//    report_abel_inversion_centers ();
//    report_abel_axes ();
//    report_abel_improper_axes ();
//    report_abel_planes (); /* order as in report_symmetry_elements_verbose() */
//  }
//  report_abel_symmetry_elements_brief ();
//  identify_point_group ();
//  symmetrize_axes (1);		/* only abelian symmetry elements used */
//  symmetrize_improper_axes (1);
//  symmetrize_inversion_centers (1);
//  symmetrize_planes (1);  /* best in this reversed order as planes preferred */
//  canonize_atoms (); 
//  printf ("--------------------------------------------\n");
//  printf ("Abelian symmetrized coordinates of all atoms\n");
//  printf ("--------------------------------------------\n");
//  report_all_atoms ();
//  mark_symmetric_planes (1);		/* particular order not important */
//  mark_symmetric_axes (1);
//  mark_symmetric_inversion_centers (1);
//  mark_symmetric_improper_axes (1);
//  printf ("--------------------------------------------\n");
//  printf ("Coordinates of abelian symmetry unique atoms\n");
//  printf ("--------------------------------------------\n");
//  report_all_atoms ();
//}
//
///*************************************************************************
// * End SGE code
// *************************************************************************/
//
//exit( EXIT_SUCCESS ) ;
//}
