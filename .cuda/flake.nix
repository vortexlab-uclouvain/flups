{
  description = "FLUPS dev shell for Nvidia GPUs";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  };

  outputs =
    { self, nixpkgs }:
    let
      system = "x86_64-linux";

      # Packages with CUDA support enabled
      pkgs = import nixpkgs {
        inherit system;
        config = {
          allowUnfree = true;
          cudaSupport = true;
          cudaCapabilities = [ "8.6" ];

          packageOverrides = pkgs: {
            # Pin CUDA version to 13.0 (at the time of writing, cudaPacakges uses 12.8)
            cudaPackages = pkgs.cudaPackages_13;
            # For some reason CUDA support doesn't work if this flag is set...
            mpi = pkgs.mpi.overrideAttrs (old: {
              configureFlags = pkgs.lib.filter (
                flag: builtins.match "--with-cuda-libdir.*" flag == null
              ) old.configureFlags or [ ];
            });
          };
        };
      };
    in
    {
      packages = {
        h3lpr = pkgs.stdenv.mkDerivation {
          name = "h3lpr";
          src = fetchGit {
            url = "git@github.com:van-Rees-Lab/h3lpr.git";
            rev = "59c63682a1eba1c20de06ae2c6dbfb8aa16610e2";
          };
          buildInputs = [ pkgs.mpi ];
          buildPhase = ''
            echo "CXX = mpic++" > make_arch/make.temp
            echo "CC = mpicc" >> make_arch/make.temp
            echo "CXXFLAGS := -O3 -std=c++17 -DCOLOR_PROF " >> make_arch/make.temp
            echo "CCFLAGS := -O3 -std=c99 -DCOLOR_PROF" >> make_arch/make.temp
            echo "LDFLAGS += -lstdc++ -lm" >> make_arch/make.temp
            ARCH_FILE=make_arch/make.temp make -j4
          '';
          installPhase = ''
            ARCH_FILE=make_arch/make.temp PREFIX=$out make install -j4
          '';
        };
      };

      devShells.${system}.default = pkgs.mkShell rec {
        buildInputs = with pkgs; [
          # CUDA stuff
          cudaPackages.cudatoolkit
          cudaPackages.nsight_systems
          cudaPackages.nsight_compute
          cudaPackages.cuda_gdb
          # Code quality tools
          pkgs.pre-commit
          cpplint
          statix
          nixfmt
          # Dev tools
          gdb
          # FLUPS deps
          mpi
          hdf5-mpi
          self.packages.h3lpr
        ];

        shellHook = ''
          export CUDA_PATH=${pkgs.cudatoolkit}
          export LD_LIBRARY_PATH=${
            pkgs.lib.makeLibraryPath (
              [
                "/run/opengl-driver" # Needed to find libGL.so
              ]
              ++ buildInputs
            )
          }:$LD_LIBRARY_PATH

          # Set environment variables for Murphy compilation
          export ARCH_FILE=make_arch/make.nix_cuda
          export CUDA_DIR=${pkgs.cudaPackages.cudatoolkit}
          export MPI_DIR=${pkgs.mpi}
          export H3LPR_DIR=${self.packages.h3lpr}
          export HDF5_DIR=${pkgs.hdf5-mpi}

          export OMPI_CC=nvcc
          export OMPI_CXX=nvcc

          # Print some info
          echo "FLUPS dev shell for Nvidia GPUs loaded"
          echo "ARCH_FILE set to make_arch/make.nix_cuda"
        '';
      };
    };
}
