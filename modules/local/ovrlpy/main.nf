/*
 * OVRLPY: Vertical signal integrity of a Xenium bundle's transcripts.
 *
 * Low integrity marks places where cells overlap in z, so no 2D segmentation
 * can separate their transcripts. The maps depend only on the transcripts;
 * cell_integrity.csv summarises them per cell of the bundle's segmentation.
 *
 * Container: params.extensions_container (an image with ovrlpy).
 *
 * Input:
 *   - meta:   Sample metadata map
 *   - bundle: Xenium bundle (the raw one, or one made by import-segmentation)
 *
 * Output:
 *   - results: folder with signal_integrity.png, integrity_map.npy,
 *              signal_map.npy, grid.json, doublets.csv, cell_integrity.csv
 */
process OVRLPY {
    tag "${meta.id}"
    label 'process_high'

    // container: params.extensions_container (conf/modules.config)

    input:
    tuple val(meta), path(bundle, stageAs: 'bundle')

    output:
    tuple val(meta), path("${prefix}"), emit: results
    tuple val("${task.process}"), val('ovrlpy'), eval("python3 -c 'import importlib.metadata as m; print(m.version(\"ovrlpy\"))'"), topic: versions, emit: versions_ovrlpy

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    ovrlpy_integrity.py \\
        --bundle bundle \\
        --outdir ${prefix} \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/signal_integrity.png ${prefix}/cell_integrity.csv
    """
}
