#if (MULT==0)
void Solver::dothemagic_std_complex_p1(double *data) {
#elif (MULT==1)
void Solver::dothemagic_std_complex_m1(double *data) {
#elif (MULT==2)
void Solver::dothemagic_std_complex_pi(double *data) {
#elif (MULT==3)
void Solver::dothemagic_std_complex_mi(double *data) {
#endif

    BEGIN_FUNC;
    int cdim = _ndim-1; // get current dim
    FLUPS_CHECK(_topo_hat[cdim]->nf() == 2, "The topo_hat[2] (field) has to be complex", LOCATION);
    // get the axis
    const int ax0 = _topo_hat[cdim]->axis();
    const int ax1 = (ax0 + 1) % 3;
    const int ax2 = (ax0 + 2) % 3;
    // get the factors
    const double         normfact = _normfact;
    opt_double_ptr       mydata   = data;
    const opt_double_ptr mygreen  = _green;
    FLUPS_ASSUME_ALIGNED(mydata,FLUPS_ALIGNMENT);
    FLUPS_ASSUME_ALIGNED(mygreen,FLUPS_ALIGNMENT);
    {
        const size_t onmax_f = _topo_hat[cdim]->nloc(ax1) * _topo_hat[cdim]->nloc(ax2) * _topo_hat[cdim]->lda();
        const size_t onmax_g = _topo_hat[cdim]->nloc(ax1) * _topo_hat[cdim]->nloc(ax2);
        const size_t inmax   = _topo_hat[cdim]->nloc(ax0);
        const int    nmem[3] = {_topo_hat[cdim]->nmem(0), _topo_hat[cdim]->nmem(1), _topo_hat[cdim]->nmem(2)};

        FLUPS_CHECK(FLUPS_ISALIGNED(mygreen) && (nmem[ax0] * _topo_hat[cdim]->nf() * sizeof(double)) % FLUPS_ALIGNMENT == 0, "please use FLUPS_ALIGNMENT to align the memory", LOCATION);
        FLUPS_CHECK(FLUPS_ISALIGNED(mydata) && (nmem[ax0] * _topo_hat[cdim]->nf() * sizeof(double)) % FLUPS_ALIGNMENT == 0, "please use FLUPS_ALIGNMENT to align the memory", LOCATION);

        // do the loop
#pragma omp parallel for default(none) proc_bind(close) schedule(static) firstprivate(onmax_f,onmax_g, inmax, nmem, mydata, mygreen, normfact, ax0)
        for (int io = 0; io < onmax_f; io++) {
            opt_double_ptr greenloc = mygreen + collapsedIndex(ax0, 0, io%onmax_g, nmem, 2); //lda of Green is only 1
            opt_double_ptr dataloc  = mydata + collapsedIndex(ax0, 0, io, nmem, 2);
            FLUPS_ASSUME_ALIGNED(dataloc,FLUPS_ALIGNMENT);
            FLUPS_ASSUME_ALIGNED(greenloc,FLUPS_ALIGNMENT);
            for (size_t ii = 0; ii < inmax; ii++) {
                const double a = dataloc[ii * 2 + 0];
                const double b = dataloc[ii * 2 + 1];
                const double c = greenloc[ii * 2 + 0];
                const double d = greenloc[ii * 2 + 1];
                // update the values
#if (MULT == 0)
                dataloc[ii * 2 + 0] = normfact * (a * c - b * d);
                dataloc[ii * 2 + 1] = normfact * (a * d + b * c);
#elif (MULT == 1)
                dataloc[ii * 2 + 0] = -normfact * (a * c - b * d);
                dataloc[ii * 2 + 1] = -normfact * (a * d + b * c);
#elif (MULT == 2)
                dataloc[ii * 2 + 0] = -normfact * (a * d + b * c);
                dataloc[ii * 2 + 1] = normfact * (a * c - b * d);
#elif (MULT == 3)
                dataloc[ii * 2 + 0] = normfact * (a * d + b * c);
                dataloc[ii * 2 + 1] = -normfact * (a * c - b * d);
#endif
            }
        }
    }
    END_FUNC;
}
