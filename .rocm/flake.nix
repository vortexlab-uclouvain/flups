{
  description = "FLUPS dev shell for AMD GPUs";

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
          rocmSupport = true;
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
          buildInputs = [ pkgs.rocmPackages.mpi ];
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
          # ROCm stuff
          rocmPackages.mpi
          rocmPackages.clr
          rocmPackages.clang
          rocmPackages.hipcc
          rocmPackages.rocfft
          rocmPackages.rocblas
          rocmPackages.rocgdb
          rocmPackages.rocm-core
          rocmPackages.rocm-device-libs
          # Code quality tools
          pkgs.pre-commit
          cpplint
          statix
          nixfmt
          # FLUPS deps
          hdf5-mpi
          self.packages.h3lpr
        ];

        shellHook = ''
          export LD_LIBRARY_PATH=${
            pkgs.lib.makeLibraryPath (
              [
                "/run/opengl-driver"
              ]
              ++ buildInputs
            )
          }:$LD_LIBRARY_PATH

          # Set environment variables for Murphy compilation
          export ARCH_FILE=make_arch/make.nix_rocm
          export MPI_DIR=${pkgs.rocmPackages.mpi}
          export H3LPR_DIR=${self.packages.h3lpr}
          export HDF5_DIR=${pkgs.hdf5-mpi}
          export ROCM_DIR=${pkgs.rocmPackages.rocm-core}
          export ROCFFT_DIR=${pkgs.rocmPackages.rocfft}
          export ROCBLAS_DIR=${pkgs.rocmPackages.rocblas}

          # Print some info
          echo "FLPUS dev shell for AMD GPUs loaded"
          echo "ARCH_FILE set to make_arch/make.nix_rocm"
        '';
      };
    };
}
