/*
 * INSTANSEG: Nuclei and cell masks from the Xenium morphology image with InstanSeg.
 *
 * Reads every channel in the bundle's morphology_focus/ (or the samplesheet's
 * image), segments one pyramid level (params.instanseg_level; 1 = half
 * resolution, close to the 0.5 um pixels InstanSeg was trained at) and writes
 * uint32 label masks at full morphology resolution, in pixel units: the input
 * Xenium Ranger import-segmentation expects.
 *
 * Runs on a GPU when params.use_gpu is true (label process_gpu_single).
 * Container: params.extensions_container, an image with instanseg-torch and
 * its model weights built in (compute nodes cannot download them).
 *
 * Input:
 *   - meta:   Sample metadata map
 *   - bundle: Xenium output bundle
 *   - image:  Morphology image to use instead of the bundle's ([] for none)
 *
 * Output:
 *   - nuclei: uint32 nuclei mask TIFF
 *   - cells:  uint32 cell mask TIFF
 */
process INSTANSEG {
    tag "${meta.id}"
    label 'process_gpu_single'

    // container: params.extensions_container (conf/modules.config)

    input:
    tuple val(meta), path(bundle, stageAs: 'bundle'), path(image)

    output:
    tuple val(meta), path("${prefix}_nuclei.tif"), emit: nuclei
    tuple val(meta), path("${prefix}_cells.tif"),  emit: cells
    tuple val("${task.process}"), val('instanseg'), eval("python3 -c 'import importlib.metadata as m; print(m.version(\"instanseg-torch\"))'"), topic: versions, emit: versions_instanseg

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def image_arg = image ? "--image ${image}" : ''
    """
    instanseg_segment.py \\
        --bundle bundle \\
        ${image_arg} \\
        --nuclei ${prefix}_nuclei.tif \\
        --cells ${prefix}_cells.tif \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_nuclei.tif ${prefix}_cells.tif
    """
}
