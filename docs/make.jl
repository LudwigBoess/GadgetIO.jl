using Documenter
using GadgetIO

makedocs(
    sitename = "GadgetIO.jl",
    format = Documenter.HTML(),
    modules = [GadgetIO],
    pages = [
            "Table of Contents" => "index.md",
            "Install"           => "install.md",
            "File Infos"        => "file_infos.md",
            "Reading Snapshots" => "read_snapshots.md",
            "Reading Subfind"   => "read_subfind.md",
            "Writing Data"      => "write_data.md",
            "API reference"     => "api.md"
            ]
        )

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/LudwigBoess/GadgetIO.jl.git"
)
