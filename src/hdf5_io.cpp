/**
 * @file hdf5_io.cpp
 * @author Denis-Gabriel Caprace and Thomas Gillis
 * @brief 
 * @version
 * @date 2019-08-20
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 * 
 */

#include "hdf5_io.hpp"

void hdf5_dump(const Topology *topo, const string filename, const double *data) {
#ifdef DUMP_H5
    xmf_write(topo, filename, "data");
    hdf5_write(topo, filename, "data", data);
#endif
}

/**
 * @brief writes the hdf5 file referenced in an xmf file
 * 
 * If the topology is real it creates one dataset `attribute`.
 * If the topology is complex it creates two datasets `attribute_real` and `attribute_imag`
 * 
 * 
 * @param topo topology of data being exported
 * @param filename the filename for export without its extension
 * @param attribute the name of the attribute to be added
 * @param data the array containing the data associated to the #topo
 * 
 * ---------------------------------------
 * Inspired from https://portal.hdfgroup.org/display/HDF5/Writing+by+Chunk+in+PHDF5
 * 
 * **General comments:**
 * A dataset is build on a dataspace (describing the datalayout) and a datatype (describing one element)
 * 
 * @warning
 * - the fastest rotating index is located at the last positions in the hdf5 C convention: `(axis+2 , axis+1 , axis)`
 * - the stride in the hyperslabs functions is the stride of 2 successive blocks!!
 * 
 * **Compression:**
 * The compression used is GNU gzip. The compression is characterize by the level of compression:
 * 
 * - `0`: No compression
 * - `1`: Best compression speed, least compression ratio
 * - `2` to `8`: compression speed degrades, compression ratio improves
 * - `9`: Slowest compression speed, best compression ratio
 * 
 * 
 * @todo for the moment no compression is set since it fails when changing the type from double to float
 * 
 * We do the following steps:
 * 
 */
void hdf5_write(const Topology *topo, const string filename, const string attribute, const double *data) {
    BEGIN_FUNC

    int mpi_size, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    hid_t  file_id;                         // file id
    hid_t  filespace_real, filespace_imag;  //dataspaces
    hid_t  fileset_real, fileset_imag;      // datasets
    hid_t  memspace;
    hid_t  plist_id; /* property list identifier */
    herr_t status;   // error code

    string extFilename = "data/" + filename + ".h5";

    // compression level
    // int complevel = 0;

    const int ax0 = topo->axis();
    const int ax1 = (ax0 + 1) % 3;
    const int ax2 = (ax0 + 2) % 3;

    //-------------------------------------------------------------------------
    /** - Create a new file collectively  */
    //-------------------------------------------------------------------------
    // setup the property list for file access (property list = option list)
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    // do some magic
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    // create the file ID
    file_id = H5Fcreate(extFilename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    if (file_id < 0) FLUPS_ERROR("Failed to open the file.");
    // close the property list
    H5Pclose(plist_id);

    //-------------------------------------------------------------------------
    /** - Create the memory dataspace  */
    //-------------------------------------------------------------------------
    // dataspace
    hsize_t memsize[3] = {(hsize_t)topo->nloc(ax2), (hsize_t)topo->nloc(ax1), (hsize_t)(topo->nloc(ax0) * topo->nf())};
    memspace           = H5Screate_simple(3, memsize, NULL);

    //-------------------------------------------------------------------------
    /** - Create the file dataspace and dataset  */
    //-------------------------------------------------------------------------
    // the file information is given by the global size
    hsize_t field_dims[3] = {(hsize_t)topo->nglob(ax2), (hsize_t)topo->nglob(ax1), (hsize_t)topo->nglob(ax0)};  // field global dimensions

    // setup the property list = option list
    plist_id           = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chk_dim[3] = {8, 8, 8};
    H5Pset_chunk(plist_id, 3, chk_dim);

    // create dataset and dataspace
    if (!topo->isComplex()) {
        filespace_real = H5Screate_simple(3, field_dims, NULL);
        fileset_real   = H5Dcreate(file_id, attribute.c_str(), H5T_NATIVE_FLOAT, filespace_real, H5P_DEFAULT, plist_id, H5P_DEFAULT);
        H5Sclose(filespace_real);
    } else {
        string realname = attribute + "_real";
        filespace_real  = H5Screate_simple(3, field_dims, NULL);
        fileset_real    = H5Dcreate(file_id, realname.c_str(), H5T_NATIVE_FLOAT, filespace_real, H5P_DEFAULT, plist_id, H5P_DEFAULT);

        string imagname = attribute + "_imag";
        filespace_imag  = H5Screate_simple(3, field_dims, NULL);
        fileset_imag    = H5Dcreate(file_id, imagname.c_str(), H5T_NATIVE_FLOAT, filespace_imag, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    }
    // close property list
    H5Pclose(plist_id);

    //-------------------------------------------------------------------------
    /** - select the hyperslab inside the file dataset (=writting location)  */
    //-------------------------------------------------------------------------
    // get the offset from topo
    int topo_offset[3];
    get_istart_glob(topo_offset, topo);

    // compute some memory quantities
    hsize_t count[3]  = {1, 1, 1};                                                                          // how many blocks to write
    hsize_t stride[3] = {1, 1, 1};                                                                          // distance between 2 blocks
    hsize_t block[3]  = {(hsize_t)topo->nloc(ax2), (hsize_t)topo->nloc(ax1), (hsize_t)topo->nloc(ax0)};     // the block size = the local size
    hsize_t offset[3] = {(hsize_t)topo_offset[ax2], (hsize_t)topo_offset[ax1], (hsize_t)topo_offset[ax0]};  // offset in the file

    // get the hyperslab within the dataset
    if (!topo->isComplex()) {
        filespace_real = H5Dget_space(fileset_real);
        status         = H5Sselect_hyperslab(filespace_real, H5S_SELECT_SET, offset, stride, count, block);
        FLUPS_CHECK0(status >= 0, "Failed to select hyperslab in dataset.")
    } else {
        filespace_real = H5Dget_space(fileset_real);
        status         = H5Sselect_hyperslab(filespace_real, H5S_SELECT_SET, offset, stride, count, block);
        FLUPS_CHECK0(status >= 0, "Failed to select real hyperslab in dataset.")
        filespace_imag = H5Dget_space(fileset_imag);
        status         = H5Sselect_hyperslab(filespace_imag, H5S_SELECT_SET, offset, stride, count, block);
        FLUPS_CHECK0(status >= 0, "Failed to select complex hyperslab in dataset.")
    }

    //-------------------------------------------------------------------------
    /** - do the writting  */
    //-------------------------------------------------------------------------
    // set the property list
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    // get the data counts
    hsize_t memblock[3] = {1, 1, 1};
    hsize_t memcount[3] = {(hsize_t)topo->nloc(ax2), (hsize_t)topo->nloc(ax1), (hsize_t)topo->nloc(ax0)};

    if (!topo->isComplex()) {
        // get the hyperslab within the dataset
        hsize_t memoffset[3] = {0, 0, 0};  // offset in memory
        hsize_t memstride[3] = {1, 1, 1};
        status               = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memoffset, memstride, memcount, memblock);
        FLUPS_CHECK0(status >= 0, "Failed to select hyperslab in memmory.")
        status = H5Dwrite(fileset_real, H5T_NATIVE_DOUBLE, memspace, filespace_real, plist_id, data);
        FLUPS_CHECK0(status >= 0, "Failed to write hyperslab to file.")
    }

    if (topo->isComplex()) {
        // stride is 2 for complex numbers
        hsize_t memstride[3] = {1, 1, 2};
        // real part
        hsize_t memoffset[3] = {0, 0, 0};
        status               = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memoffset, memstride, memcount, memblock);
        FLUPS_CHECK0(status >= 0, "Failed to select real hyperslab in memmory.")
        status = H5Dwrite(fileset_real, H5T_NATIVE_DOUBLE, memspace, filespace_real, plist_id, data);
        FLUPS_CHECK0(status >= 0, "Failed to write real part hyperslab to file.")

        // imaginary part
        memoffset[2] = 1;  // set an offset on the fastest rotating index
        status       = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memoffset, memstride, memcount, memblock);
        FLUPS_CHECK0(status >= 0, "Failed to select imag hyperslab in memmory.")
        status = H5Dwrite(fileset_imag, H5T_NATIVE_DOUBLE, memspace, filespace_imag, plist_id, data);
        FLUPS_CHECK0(status >= 0, "Failed to write imaginary part hyperslab to file.")
    }

    //-------------------------------------------------------------------------
    /** - close everything */
    //-------------------------------------------------------------------------
    H5Dclose(fileset_real);
    H5Sclose(filespace_real);
    H5Sclose(memspace);
    if (topo->isComplex()) {
        H5Dclose(fileset_imag);
        H5Sclose(filespace_imag);
    }

    H5Pclose(plist_id);
    H5Fclose(file_id);

    return;
}

/**
 * @brief writes a xmf file, readable by a XDMF viewer
 * 
 * @param topo topology of data being exported
 * @param filename the filename for export without its extension
 * @param attribute the name of the attribute to be added
 * 
 * -----------------------------
 * inspired from http://visitusers.org/index.php?title=Using_XDMF_to_read_HDF5
 * see also http://www.xdmf.org/index.php/XDMF_Model_and_Format
 */
void xmf_write(const Topology *topo, const string filename, const string attribute) {
    BEGIN_FUNC

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const int ax0 = topo->axis();
    const int ax1 = (ax0 + 1) % 3;
    const int ax2 = (ax0 + 2) % 3;

    FILE * xmf         = 0;
    string extFilename = "data/" + filename + ".xmf";

    if (rank == 0) {
        xmf = fopen(extFilename.c_str(), "w");

        fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
        fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
        fprintf(xmf, "<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"3.0\">\n");
        fprintf(xmf, "<!-- Kindly generated by poisson -->\n");
        fprintf(xmf, " <Domain>\n");
        fprintf(xmf, "   <Grid GridType=\"Uniform\">\n");

        // the size of the grid has to be +1 since the 3DCoRectMesh is defined as vertex-centered
        fprintf(xmf, "     <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"%d %d %d\"/>\n", topo->nglob(ax2) + 1, topo->nglob(ax1) + 1, topo->nglob(ax0) + 1);

        fprintf(xmf, "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n");
        fprintf(xmf, "           <DataItem Name=\"Origin\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
        fprintf(xmf, "              %10.8f %10.8f %10.8f\n", 0.0, 0.0, 0.0);
        fprintf(xmf, "           </DataItem>\n");
        fprintf(xmf, "           <DataItem Name=\"Spacing\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
        fprintf(xmf, "              %10.8f %10.8f %10.8f\n", 1.0, 1.0, 1.0);
        fprintf(xmf, "           </DataItem>\n");
        fprintf(xmf, "     </Geometry>\n");

        if (!topo->isComplex()) {
            fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n", attribute.c_str());
            fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", topo->nglob(ax2), topo->nglob(ax1), topo->nglob(ax0));
            fprintf(xmf, "        %s.h5:/%s\n", filename.c_str(), attribute.c_str());  //<-------------------------------
            fprintf(xmf, "       </DataItem>\n");
            fprintf(xmf, "     </Attribute>\n");
        } else {
            // real part
            string realname = attribute + "_real";
            fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n", realname.c_str());
            fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", topo->nglob(ax2), topo->nglob(ax1), topo->nglob(ax0));
            fprintf(xmf, "        %s.h5:/%s\n", filename.c_str(), realname.c_str());  //<-------------------------------
            fprintf(xmf, "       </DataItem>\n");
            fprintf(xmf, "     </Attribute>\n");
            // imag part
            string imagname = attribute + "_imag";
            fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n", imagname.c_str());
            fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", topo->nglob(ax2), topo->nglob(ax1), topo->nglob(ax0));
            fprintf(xmf, "        %s.h5:/%s\n", filename.c_str(), imagname.c_str());  //<-------------------------------
            fprintf(xmf, "       </DataItem>\n");
            fprintf(xmf, "     </Attribute>\n");
        }

        //For Node centered:
        // fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n", scalarName);
        // fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", NY, NX);
        // fprintf(xmf, "        %s.h5:/data\n", filename); //<-------------------------------
        // fprintf(xmf, "       </DataItem>\n");
        // fprintf(xmf, "     </Attribute>\n");

        //For Vector:
        // fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Vector\" Center=\"Node\">\n", scalarName);
        // fprintf(xmf, "       <DataItem Dimensions=\"%d %d  3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", NY, NX);
        // fprintf(xmf, "        xdmf2d.h5:/data\n"); //<-------------------------------
        // fprintf(xmf, "       </DataItem>\n");
        // fprintf(xmf, "     </Attribute>\n");
        //Also possible: Tensor (3x3), Tensor6 or even Matrix of NxM

        fprintf(xmf, "   </Grid>\n");
        fprintf(xmf, " </Domain>\n");
        fprintf(xmf, "</Xdmf>\n");
        fclose(xmf);
    }
}

/**
 * @brief test the hdf5 dumping for real and complex numbers
 * 
 */
void hdf5_dumptest() {
    BEGIN_FUNC

    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    const int nglob[3] = {8, 8, 8};
    const int nproc[3] = {comm_size, 1, 1};

    //===========================================================================
    // real numbers
    Topology *topo = new Topology(0, nglob, nproc, false);

    double *data = (double *)fftw_malloc(sizeof(double *) * topo->locmemsize());

    for (int i2 = 0; i2 < topo->nloc(2); i2++) {
        for (int i1 = 0; i1 < topo->nloc(1); i1++) {
            for (int i0 = 0; i0 < topo->nloc(0); i0++) {
                size_t id    = localindex_xyz(i0, i1, i2, topo);
                data[id + 0] = id;
            }
        }
    }
    // try the dump
    hdf5_dump(topo, "test_real", data);

    fftw_free(data);
    delete (topo);

    //===========================================================================
    // create a real topology
    topo = new Topology(0, nglob, nproc, true);

    data = (double *)fftw_malloc(sizeof(double *) * topo->locmemsize());

    for (int i2 = 0; i2 < topo->nloc(2); i2++) {
        for (int i1 = 0; i1 < topo->nloc(1); i1++) {
            for (int i0 = 0; i0 < topo->nloc(0); i0++) {
                size_t id    = localindex_xyz(i0, i1, i2, topo);
                data[id + 0] = (double)id;
                data[id + 1] = 0.0;  //-data[id + 0];
            }
        }
    }
    // try the dump
    hdf5_dump(topo, "test_complex", data);
    fftw_free(data);
    delete (topo);
}