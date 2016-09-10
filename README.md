# Yet Another Pipeline System (YAPS) 2 _(COSMOS2-based)_

_A work in progress.  Not really ready for prime time.  This is a pure idea experiment and code may be abandoned in the future._

_[yaps][0] is based on [ruffus][1].  This repository, [yaps2][3], is based on a tweaked [COSMOS2][2]._

## Example execution

    yaps2 postvqsr --workspace=/path/to/workspace/name/test-pipeline --input-vcfs=inputs.txt --project-name="test-pipeline" --timeout=300 --log-level=INFO

## Deployment

    pip install --process-dependency-links git+https://github.com/indraniel/yaps2.git
    # pygraphviz needs some extra install help
    pip install pygraphviz --install-option="--include-path=/usr/include/graphviz" --install-option="--library-path=/usr/lib/graphviz/" --force-reinstall --upgrade

## Development

    git clone https://github.com/indraniel/yaps2.git
    cd yaps2
    virtualenv venv
    source venv/bin/activate
    pip install --process-dependency-links -e .
    # pygraphviz needs some extra install help
    pip install pygraphviz --install-option="--include-path=/usr/include/graphviz" --install-option="--library-path=/usr/lib/graphviz/" --force-reinstall --upgrade
    # < do development work >
     
    # test
    yaps2 postvqsr --workspace=/path/to/workspace/name/test-pipeline --input-vcfs=inputs.txt --project-name="test-pipeline" --timeout=300 --log-level=INFO

    # clean up dev workspace
    make clean

## Misc. Notes

* Currently only the `postvqsr` pipeline exists.  See `yaps2 --help` and/or `yaps postvqsr --help` for more information on the available commands/pipelines and options.
* `--input-vcfs` is a file containing a tab-separated list of `*.vcf.gz` files in `<CHROM>\t<VCF.GZ FILE>` format

[0]: https://github.com/indraniel/yaps
[1]: http://www.ruffus.org.uk/
[2]: https://github.com/indraniel/COSMOS2/tree/enable-lsf-rebase
[3]: https://github.com/indraniel/yaps2
