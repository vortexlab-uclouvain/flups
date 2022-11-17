/**
 * @file FFTW_plan_dim.hpp
 * @copyright Copyright © Université catholique de Louvain (UCLouvain), Belgique 
 *      See LICENSE file in top-level directory
*/

#ifndef FFTW_PLAN_DIM_HPP
#define FFTW_PLAN_DIM_HPP

#include <array>
#include <tuple>

#include "Topology.hpp"
#include "defines.hpp"
#include "fftw3.h"

/**
 * @brief A FFTW plan in one dimension
 *
 */
class FFTW_plan_dim {
   public:
    /**
     * @brief PlanType is the type of plan considered and is computed as the sum of both BoundaryType variables
     *
     * The integer value associated gives is the priority of processing.
     * We first have to do the real to real transforms, then the padded real to real (mix direction = unbounded + boundary condition),
     * then the periodic (DFT) directions and finally the padded periodic boundary condition.
     * This order is chosen in order to reduce the computational cost.
     *
     * If a multi-dimension FFT is asked, one plan is created for each dimension as it may be different.
     * If the plans are the same, we keep the first plan issued
     */
    enum PlanType {
        SYMSYM = 2, /**< type real 2 real (DCT / DST) : EE (0) , EO/OE (1) , OO (2) */
        MIXUNB = 5, /**< type unbounded and a symetry condition: UE/EU (4) , UO/OU (5) */
        PERPER = 6, /**< type periodic - periodic: PERPER (6) */
        UNBUNB = 8, /**< type fully unbounded UU (8) */
        EMPTY  = 18 /**< type empty, i.e. this direction is not used */
    };

    /**
     * @brief Determines which post processing operation must be performed on the output of the fft forward/backward
     *
     * the corrections combines 3 basic operations and therefore are encrypted on 3 bytes:
     * - byte 0 obtained as `(correction>>0)%2` overwrites the first data of the transform
     * - byte 1 obtained as `(correction>>1)%2` overwrites the last data of the transform
     * - byte 2 obtained as `(correction>>1)%2` copy the value of the first data in the last data
     *
     * The actual correction of a planned is asigned using the sum operator:
     *       correction = NULL_LAST_POINT = 2 
     * will lead to
     * - (2>>0)%2 = 0 -> NO first point correction
     * - (2>>1)%2 = 1 -> YES last point correction
     * - (2>>2)%2 = 1 -> NO periodic correction
     */
    enum PlanPostproType {
        POSTPRO_NONE     = 0,  // no corrections
        NULL_FIRST_POINT = 1,  // byte 0, obtained as = 1<<0
        NULL_LAST_POINT  = 2,  // byte 1, obtained as = 1<<1
        ENFORCE_PERIOD   = 4,  // byte 2, obtained as = 2<<1
    };

   protected:
    const int  lda_;                /**<@brief the dimension of the solver */
    const bool isGreen_;            /**< @brief boolean is true if this plan is for a Green's function */
    const int  dimID_;              /**< @brief the dimension of the plan in the field reference */
    const int  sign_;               /**< @brief FFT_FORWARD (-1) or FFT_BACKWARD(+1) */

    bool   isr2c_       = false;    /**< @brief is this plan the one that changes to complex?*/
    bool   isSpectral_  = false;    /**< @brief indicate if the Green's function has to be done spectrally (leading to a helmolz problem) */
    int    fftw_stride_ = 0;        /**<@brief the memory space between two ffts */
    int    howmany_     = 0;        /**<@brief the number of FFT's to do */
    
    int    n_out_       = 1;        /**< @brief the number of element coming out of the transform. When dealing with vector, this number must be constant throughout the different components*/
    int    fieldstart_  = 0;        /**< @brief the starting index for the field copy in the direction of the plan*/
    double symstart_    = 0.0;      /**< @brief the first index to be copied for the symmetry done on the Green's function, set to 0 if no symmetry is needed*/
    double volfact_     = 1.0;      /**< @brief volume factor*/
    double normfact_    = 1.0;      /**< @brief factor you need to multiply to get the transform on the right scaling*/
    double kfact_       = 0.0;      /**< @brief multiplication factor to have the correct k numbers*/
    double koffset_     = 0.0;      /**< @brief additive factor to have the correct k numbers*/

    int* n_in_          = NULL;     /**< @brief the number of element in the transform, i. e. given to fftw calls*/
    int* fftwstart_in_  = NULL;     /**< @brief the starting index for the input field to be given to FFTW functions*/
    int* fftwstart_out_ = NULL;     /**< @brief the starting index for the output field to be given to FFTW functions*/

    PlanType       type_;                        /**< @brief type of this plan, see #PlanType*/
    BoundaryType*  bc_[2]        = {NULL, NULL}; /**< @brief boundary condition for the ith component [0][i]=LEFT/MIN - [1][i]=RIGHT/MAX*/
    int*           postpro_type_ = NULL;         /**< @brief correction type of this plan, see #PlanPostproType*/
    bool*          imult_        = NULL;         /**< @brief boolean indicating that we have to multiply by (-i) in forward and (i) in backward*/
    fftw_r2r_kind* kind_         = NULL;         /**< @brief kind of transfrom to perform (used by r2r and mix plan only)*/
    fftw_plan*     plan_         = NULL;         /**< @brief the array of FFTW plan*/

   public:
    FFTW_plan_dim(const int lda, const int dimID, const double h[3], const double L[3], BoundaryType* mybc[2], const int sign, const bool isGreen);
    virtual ~FFTW_plan_dim();  // Virtual, following http://www.gotw.ca/publications/mill18.htm

    void init(const int size[3], const bool isComplex);

    void allocate_plan(const Topology* topo, double* data);
    void execute_plan(const Topology* topo, double* data) const;
    void postprocess_plan(const Topology*, double* data);

    /**
     * @name Getters - return the value
     *
     */
    /**@{ */
    inline bool   isSpectral() const { return isSpectral_; }
    inline bool   isr2c() const { return isr2c_; }
    inline bool   imult(const int lia) const { return imult_[lia]; }
    inline bool   isr2c_doneByFFT() const { return isr2c_ && (!isSpectral_); }
    inline int    dimID() const { return dimID_; }
    inline int    type() const { return type_; }
    inline double symstart() const { return symstart_; }
    inline double normfact() const { return normfact_; }
    inline double volfact() const { return volfact_; }
    inline double kfact() const { return kfact_; }
    inline double koffset() const { return koffset_; }
    inline void   get_outsize(int* size) const { size[dimID_] = n_out_; };
    inline void   get_fieldstart(int* start) const { start[dimID_] = fieldstart_; };
    inline void   get_isNowComplex(bool* isComplex) const { (*isComplex) = (*isComplex) || isr2c_; };
    /**@} */

    /**
     * @name correction helper functions
     *
     * return true if the correction required by the plan contains the associated operation
     *
     */
    /**@{ */
    bool do_reset_first_point(const int value) const { return ((value >> 0) % 2); };
    bool do_reset_last_point(const int value) const { return ((value >> 1) % 2); };
    bool do_enforce_period(const int value) const { return ((value >> 2) % 2); };
    /**@} */

    void disp();

   protected:
    void check_dataAlign_(const Topology* topo, double* data) const;

    /**
     * @name Plan allocation
     */
    /**@{ */
    void allocate_plan_real_(const Topology* topo, double* data);
    void allocate_plan_complex_(const Topology* topo, double* data);
    /**@} */

    /**
     * @name Initialization
     */
    /**@{ */
    virtual void        init_real2real_(const int size[3], bool isComplex)    = 0;
    virtual void        init_mixunbounded_(const int size[3], bool isComplex) = 0;
    virtual void        init_periodic_(const int size[3], bool isComplex)     = 0;
    virtual void        init_unbounded_(const int size[3], bool isComplex)    = 0;
    virtual void        init_empty_(const int size[3], bool isComplex)        = 0;
    virtual std::string disp_data_center() const                              = 0;
    /**@} */
};

void sort_priority(std::array<std::tuple<int, int>, 3>* priority);
void sort_plans(FFTW_plan_dim* plan[3]);
int  bc_to_types(const BoundaryType* bc[2]);

#endif