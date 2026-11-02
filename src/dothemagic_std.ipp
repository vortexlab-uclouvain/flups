/**
 * @file dothemagic_std.ipp
 * @author Thomas Gillis and Denis-Gabriel Caprace
 * @brief contains the domagic functions
 * @version
 *
 * @copyright Copyright © UCLouvain 2020
 * FLUPS is a Fourier-based Library of Unbounded Poisson Solvers.
 *
 * Copyright <2020> <Université catholique de Louvain (UCLouvain), Belgique>
 *
 * List of the contributors to the development of FLUPS, Description and complete License: see LICENSE and NOTICE files.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#if (KIND == 0)
/**
 * @brief perform the convolution for real to real cases
 *
 */
void Solver::dothemagic_std_real(double *data, stream_t stream) {
#elif (KIND == 1)
/**
 * @brief perform the convolution for complex to complex cases
 *
 */
void Solver::dothemagic_std_complex(double *data, stream_t stream) {
#endif

    BEGIN_FUNC;
    int cdim = ndim_ - 1;  // get current dim
#if (KIND == 0)
    FLUPS_CHECK(topo_hat_[cdim]->nf() == 1, "The topo_hat[2] (field) has to be complex");
#else
    FLUPS_CHECK(topo_hat_[cdim]->nf() == 2, "The topo_hat[2] (field) has to be complex");
#endif
    // get the axis
    const int nf  = topo_hat_[cdim]->nf();
    const int ax0 = topo_hat_[cdim]->axis();
    const int ax1 = (ax0 + 1) % 3;
    const int ax2 = (ax0 + 2) % 3;

    // get the norm factor
    const double normfact = normfact_;

    // get the adresses
    opt_double_ptr       mydata  = data;
    const opt_double_ptr mygreen = green_;

    // get the number of pencils for the field and green
    const size_t onmax = topo_hat_[cdim]->nloc(ax1) * topo_hat_[cdim]->nloc(ax2) * topo_hat_[cdim]->lda();
    const size_t ondim = topo_hat_[cdim]->nloc(ax1) * topo_hat_[cdim]->nloc(ax2);
    const size_t inmax = topo_hat_[cdim]->nloc(ax0);
    // get the memory details
    const size_t memdim  = topo_hat_[cdim]->memdim();
    const int    nmem[3] = {topo_hat_[cdim]->nmem(0), topo_hat_[cdim]->nmem(1), topo_hat_[cdim]->nmem(2)};
    // stride between two successive pencils (in doubles) = collapsedIndex(ax0, 0, 1, nmem, nf)
    const size_t stride_o = (size_t)nmem[ax0] * (size_t)nf;
    const int    ondim_i  = (int)ondim;

    // check the alignment
    FLUPS_CHECK(FLUPS_ISALIGNED(mygreen) && (nmem[ax0] * topo_hat_[cdim]->nf() * sizeof(double)) % FLUPS_ALIGNMENT == 0, "please use FLUPS_ALIGNMENT to align the memory");
    FLUPS_CHECK(FLUPS_ISALIGNED(mydata) && (nmem[ax0] * topo_hat_[cdim]->nf() * sizeof(double)) % FLUPS_ALIGNMENT == 0, "please use FLUPS_ALIGNMENT to align the memory");
    FLUPS_ASSUME_ALIGNED(mydata, FLUPS_ALIGNMENT);
    FLUPS_ASSUME_ALIGNED(mygreen, FLUPS_ALIGNMENT);

    // do the loop
    // one (component, pencil) per outer index, the point in the pencil is the inner (contiguous) index
    const auto op = m_backend_lambda(const int iop, const int ii) {
        // unroll the outer index
        const int io  = iop % ondim_i;
        const int lia = iop / ondim_i;

        // get the starting pointer
        opt_double_ptr greenloc = mygreen + stride_o * io;  // lda of Green is only 1
        opt_double_ptr dataloc  = mydata + lia * memdim + stride_o * io;
        FLUPS_ASSUME_ALIGNED(greenloc, FLUPS_ALIGNMENT);
        FLUPS_ASSUME_ALIGNED(dataloc, FLUPS_ALIGNMENT);

        // do the actual convolution
#if (KIND == 0)
        dataloc[ii] *= normfact * greenloc[ii];
#elif (KIND == 1)
        const double a = dataloc[ii * 2 + 0];
        const double b = dataloc[ii * 2 + 1];
        const double c = greenloc[ii * 2 + 0];
        const double d = greenloc[ii * 2 + 1];
        // update the values
        dataloc[ii * 2 + 0] = normfact * (a * c - b * d);
        dataloc[ii * 2 + 1] = normfact * (a * d + b * c);
#endif
    };
    parallel_for_2d(op, (int)onmax, (int)inmax, "dothemagic_std", stream);

    END_FUNC;
}
