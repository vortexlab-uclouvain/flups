#if (MULT==0)
void Solver::dothemagic_rot_complex_p1(double *data, double kfact[3], double koffset[3], double symstart[3], int _orderdiff) {
#elif (MULT==1)
void Solver::dothemagic_rot_complex_m1(double *data, double kfact[3], double koffset[3], double symstart[3], int _orderdiff) {
#elif (MULT==2)
void Solver::dothemagic_rot_complex_pi(double *data, double kfact[3], double koffset[3], double symstart[3], int _orderdiff) {
#elif (MULT==3)
void Solver::dothemagic_rot_complex_mi(double *data, double kfact[3], double koffset[3], double symstart[3], int _orderdiff) {
#endif

    BEGIN_FUNC;
    int cdim = _ndim-1; // get current dim
    FLUPS_CHECK(_topo_hat[cdim]->nf() == 2, "The topo_hat[2] (field) has to be complex", LOCATION);
    // get the axis
    const int ax0 = _topo_hat[cdim]->axis();
    const int ax1 = (ax0 + 1) % 3;
    const int ax2 = (ax0 + 2) % 3;
    const int nf  = _topo_hat[cdim]->nf();

    int istart[3];
    _topo_hat[cdim]->get_istart_glob(istart);

    // get the factors
    const double         normfact = _normfact;
    opt_double_ptr       mydata   = data;
    const opt_double_ptr mygreen  = _green;

    FLUPS_ASSUME_ALIGNED(mydata,FLUPS_ALIGNMENT);
    FLUPS_ASSUME_ALIGNED(mygreen,FLUPS_ALIGNMENT);
    {
        const size_t onmax   = _topo_hat[cdim]->nloc(ax1) * _topo_hat[cdim]->nloc(ax2);
        const int    nmem[3] = {_topo_hat[cdim]->nmem(0), _topo_hat[cdim]->nmem(1), _topo_hat[cdim]->nmem(2)};
        FLUPS_CHECK(FLUPS_ISALIGNED(mygreen) && (nmem[ax0] * _topo_hat[cdim]->nf() * sizeof(double)) % FLUPS_ALIGNMENT == 0, "please use FLUPS_ALIGNMENT to align the memory", LOCATION);
        FLUPS_CHECK(FLUPS_ISALIGNED(mydata) && (nmem[ax0] * _topo_hat[cdim]->nf() * sizeof(double)) % FLUPS_ALIGNMENT == 0, "please use FLUPS_ALIGNMENT to align the memory", LOCATION);

        int nloc[3] = {_topo_hat[cdim]->nloc(0), _topo_hat[cdim]->nloc(1),  _topo_hat[cdim]->nloc(2)};

        // do the loop
#pragma omp parallel for default(none) proc_bind(close) schedule(static) firstprivate(onmax, nmem, mydata, mygreen, normfact, ax0, ax1, ax2, nf, nloc, istart,symstart,kfact,koffset)
        for (int i2 = 0; i2 < nloc[ax2]; i2++) {
            for (int i1 = 0; i1 < nloc[ax1]; i1++) {
                //local indexes start
                const size_t id = localIndex(ax0, 0, i1, i2, ax0, nmem, nf, 0);

                opt_double_ptr greenloc = mygreen + collapsedIndex(ax0, 0, i1*i2, nmem, 2); //lda of Green is only 1
                opt_double_ptr dataloc  = mydata + collapsedIndex(ax0, 0, i1*i2, nmem, 2);
                FLUPS_ASSUME_ALIGNED(dataloc,FLUPS_ALIGNMENT);
                FLUPS_ASSUME_ALIGNED(greenloc,FLUPS_ALIGNMENT);
            
                for (int i0 = 0; i0 < nloc[ax0]; i0++) {
                    int is[3];
                    cmpt_symID(ax0,i0,i1,i2,istart,symstart,0,is);

                    // (symmetrized) wave number : only 1 kfact is zero
                    const double k0 = (is[ax0] + koffset[ax0]) * kfact[ax0];
                    const double k1 = (is[ax1] + koffset[ax1]) * kfact[ax1];
                    const double k2 = (is[ax2] + koffset[ax2]) * kfact[ax2];

                    const double f0r = dataloc[i0 * 2 + 0];
                    const double f0c = dataloc[i0 * 2 + 1];
                    const double f1r = dataloc[i0 * 2 + onmax + 0];
                    const double f1c = dataloc[i0 * 2 + onmax + 1];
                    const double f2r = dataloc[i0 * 2 + 2 * onmax + 0];
                    const double f2c = dataloc[i0 * 2 + 2 * onmax + 1];

                    const double rot0r = k1 * f2r - k2 * f1r;
                    const double rot0c = k1 * f2c - k2 * f1c;
                    const double rot1r = k2 * f0r - k0 * f2r;
                    const double rot1c = k2 * f0c - k0 * f2c;
                    const double rot2r = k0 * f1r - k1 * f0r;
                    const double rot2c = k0 * f1c - k1 * f0c;

                    const double gr = greenloc[i0 * 2 + 0];
                    const double gc = greenloc[i0 * 2 + 1];
 
                    // update the values
#if (MULT == 0)
                    dataloc[i0 * 2 + 0] = normfact * (rot0r * gr - rot0c * gc);
                    dataloc[i0 * 2 + 1] = normfact * (rot0r * gc + rot0c * gr);
                    dataloc[i0 * 2 + onmax + 0] = normfact * (rot1r * gr - rot1c * gc);
                    dataloc[i0 * 2 + onmax + 1] = normfact * (rot1r * gc + rot1c * gr);
                    dataloc[i0 * 2 + onmax * 2 + 0] = normfact * (rot2r * gr - rot2c * gc);
                    dataloc[i0 * 2 + onmax * 2 + 1] = normfact * (rot2r * gc + rot2c * gr);
#elif (MULT == 1)
                    dataloc[i0 * 2 + 0] = -normfact * (rot0r * gr - rot0c * gc);
                    dataloc[i0 * 2 + 1] = -normfact * (rot0r * gc + rot0c * gr);
                    dataloc[i0 * 2 + onmax + 0] = -normfact * (rot1r * gr - rot1c * gc);
                    dataloc[i0 * 2 + onmax + 1] = -normfact * (rot1r * gc + rot1c * gr);
                    dataloc[i0 * 2 + onmax * 2 + 0] = -normfact * (rot2r * gr - rot2c * gc);
                    dataloc[i0 * 2 + onmax * 2 + 1] = -normfact * (rot2r * gc + rot2c * gr);
#elif (MULT == 2)
                    dataloc[i0 * 2 + 0] = -normfact * (rot0r * gc + rot0c * gr);
                    dataloc[i0 * 2 + 1] = normfact * (rot0r * gr - rot0c * gc);
                    dataloc[i0 * 2 + onmax + 0] = -normfact * (rot1r * gc + rot1c * gr);
                    dataloc[i0 * 2 + onmax + 1] = normfact * (rot1r * gr - rot1c * gc);
                    dataloc[i0 * 2 + onmax * 2 + 0] = -normfact * (rot2r * gc + rot2c * gr);
                    dataloc[i0 * 2 + onmax * 2 + 1] = normfact * (rot2r * gr - rot2c * gc);

#elif (MULT == 3)
                    dataloc[i0 * 2 + 0] = normfact * (rot0r * gc + rot0c * gr);
                    dataloc[i0 * 2 + 1] = -normfact * (rot0r * gr - rot0c * gc);
                    dataloc[i0 * 2 + onmax + 0] = normfact * (rot1r * gc + rot1c * gr);
                    dataloc[i0 * 2 + onmax + 1] = -normfact * (rot1r * gr - rot1c * gc);
                    dataloc[i0 * 2 + onmax * 2 + 0] = normfact * (rot2r * gc + rot2c * gr);
                    dataloc[i0 * 2 + onmax * 2 + 1] = -normfact * (rot2r * gr - rot2c * gc);
#endif
                }
            }
        }
    }
    END_FUNC;
}
