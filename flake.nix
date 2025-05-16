{
  description = "Reproducible dev shell for the `NVD2` bioinformatic processing pipeline";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
      ...
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs {
          inherit system;
        };

      in
      {
        devShells.default = pkgs.mkShell {
          name = "NVD2";

          buildInputs = [
            pkgs.stdenv
            pkgs.gcc
            pkgs.curl
            pkgs.wget
            pkgs.openjdk
            pkgs.git
            pkgs.cmake
            pkgs.libxml2
            pkgs.libxslt
            pkgs.libffi
            pkgs.pixi
          ];

          shellHook = ''
            echo "🔧 Entering NVD dev shell"
            export PS1="(NVD) $PS1"
            if [ ! -d .pixi/envs/default ]; then
              echo "Pixi env not found. Running install..."
              pixi install --frozen
            fi

            export PATH="$PWD/.pixi/envs/default/bin:$PATH"
          '';
        };
      }
    );
}

