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
            "Merger Trees"      => "merger_trees.md",
            "Main Halo Finder"  => "find_main_halo.md",
            "Writing Data"      => "write_data.md",
            "Reading CPU files" => "read_cpu_files.md",
            "API reference"     => [
                "Index"              => "api.md",
                "Exported Functions" => "api_exported_functions.md",
                "Exported Types"     => "api_exported_types.md",
                "Private Functions"  => "api_private_functions.md",
                "Private Types"      => "api_private_types.md"
                ]
            ]
        )

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/LudwigBoess/GadgetIO.jl.git"
)
