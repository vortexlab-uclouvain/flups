/**
 * @file hdf5_io.cpp
 * @author Thomas Gillis
 * @brief 
 * @version
 * @date 2019-08-20
 * 
 * @copyright Copyright © UCLouvain 2019
 * 
 */

#include "hdf5_io.hpp"

/**
 * @brief writes the hdf5 file referenced in an xmf file
 * 
 * @param topo topology of data being exported
 * @param filename the filename for export without its extension
 * @param attribute the name of the attribute to be added
 * @param data the array containing the data associated to the #topo
 * 
 * ---------------------------------------
 * Inspired from https://portal.hdfgroup.org/display/HDF5/Writing+by+Chunk+in+PHDF5
 * 
 * We do the following steps:
 * 
 */
void hdf5_write(const Topology *topo, const string filename, const string attribute, const double *data)
{
    BEGIN_FUNC

    int mpi_size, mpi_rank, info;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    hid_t file_id;   // file id
    hid_t dset_id;   // dataset id
    hid_t filespace; // file dataspace identifier
    hid_t memspace;  // dataspace identifier */
    hid_t plist_id;  /* property list identifier */
    herr_t status;   // error code

    string extFilename = filename + ".h5";

    //-------------------------------------------------------------------------
    /** - Create a new file collectively  */
    //-------------------------------------------------------------------------
    // setup the property list for file access (property list = option list)
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    // do some magic
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    // create the file ID
    file_id = H5Fcreate(extFilename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    // close the property list
    H5Pclose(plist_id);

    //-------------------------------------------------------------------------
    /** - Create the dataspace (defines the organization of the elements) for the file and for the chunks  */
    //-------------------------------------------------------------------------
    // the memory information is given by local size
    // hsize_t chunk_dims[3] = {topo->nloc(0), topo->nloc(1), topo->nloc(2)};
    hsize_t chunk_dims[3] = {topo->nloc(2), topo->nloc(1), topo->nloc(0)};
    memspace = H5Screate_simple(3, chunk_dims, NULL);
    // the file information is given by the global size
    // hsize_t field_dims[3] = {topo->nglob(0), topo->nglob(1), topo->nglob(2)}; // field global dimensions
    hsize_t field_dims[3] = {topo->nglob(2), topo->nglob(1), topo->nglob(0)}; // field global dimensions
    filespace = H5Screate_simple(3, field_dims, NULL);

    //-------------------------------------------------------------------------
    /** - Create the dataset for the file and set it to the "chunked state" */
    //-------------------------------------------------------------------------
    // setup the property list = option list
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    // update the property with the chunk dimensions and the rank
    H5Pset_chunk(plist_id, 3, chunk_dims);
    // create the dataset (in floats!!)
    dset_id = H5Dcreate(file_id, attribute.c_str(), H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    // close the filespace and the property list
    H5Sclose(filespace);
    H5Pclose(plist_id);

    //-------------------------------------------------------------------------
    /** - select the hyperslab (=writting location) inside the file dataset */
    //-------------------------------------------------------------------------
    // get the offset from topo
    int topo_offset[3];
    topo->get_idstart_glob(topo_offset);

    // compute some memory quantities
    hsize_t block[3] = {topo->nloc(2), topo->nloc(1), topo->nloc(0)}; // the block size = the chunk size
    hsize_t count[3] = {1, 1, 1};                                     // howmany blocks is to write
    hsize_t stride[3] = {1, 1, 1};                                    // we take every element
    hsize_t offset[3] = {(hsize_t)topo_offset[2], (hsize_t)topo_offset[1], (hsize_t)topo_offset[0]}; // offset in memory

    // get the hyperslab within the dataset
    filespace = H5Dget_space(dset_id);
    status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, block);

    //-------------------------------------------------------------------------
    /** - write the data from the memory to the file */
    //-------------------------------------------------------------------------
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    // we write the data (double format) into the hyperslab
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);

    //-------------------------------------------------------------------------
    /** - close everything */
    //-------------------------------------------------------------------------
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
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
 */
void xmf_write(const Topology *topo, const string filename, const string attribute)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    FILE *xmf = 0;
    string extFilename = filename + ".xmf";
    
    if (rank == 0)
    {
        xmf = fopen(extFilename.c_str(), "w");

        fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
        fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
        fprintf(xmf, "<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"3.0\">\n");
        fprintf(xmf, "<!-- Kindly generated by poisson -->\n");

        fprintf(xmf, " <Domain>\n");
        fprintf(xmf, "   <Grid GridType=\"Uniform\">\n");

        // the sizxe of the grid has to be +1 since the 3DCoRectMesh is defined as vertex-centered
        fprintf(xmf, "     <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"%d %d %d\"/>\n", topo->nglob(2)+1, topo->nglob(1)+1, topo->nglob(0)+1);

        fprintf(xmf, "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n");
        fprintf(xmf, "           <DataItem Name=\"Origin\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
        fprintf(xmf, "              %10.8f %10.8f %10.8f\n", 0.0, 0.0, 0.0);
        fprintf(xmf, "           </DataItem>\n");
        fprintf(xmf, "           <DataItem Name=\"Spacing\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
        fprintf(xmf, "              %10.8f %10.8f %10.8f\n", 1.0, 1.0, 1.0);
        fprintf(xmf, "           </DataItem>\n");
        fprintf(xmf, "     </Geometry>\n");

        fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n", attribute.c_str());
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", topo->nglob(2), topo->nglob(1), topo->nglob(0));
        fprintf(xmf, "        %s.h5:/data\n", filename.c_str()); //<-------------------------------
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");

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
