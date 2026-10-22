//
// Run InstanSeg on the morphology image, then import its masks with Xenium Ranger
//

include { INSTANSEG                       } from '../../../modules/local/instanseg/main'
include { XENIUMRANGER_IMPORTSEGMENTATION } from '../../../modules/nf-core/xeniumranger/importsegmentation/main'

workflow INSTANSEG_IMPORT_SEGMENTATION {
    take:
    ch_instanseg_input         // channel: [ val(meta), ["path-to-xenium-bundle"], ["path-to-image"] or [] ]
    ch_bundle_path             // channel: [ val(meta), ["path-to-xenium-bundle"] ]
    nucleus_segmentation_only  // value: bool; import nuclei only, expanded by expansion_distance
    expansion_distance         // value: nuclear expansion distance (microns)

    main:

    ch_coordinate_space = channel.value("pixels")

    INSTANSEG(ch_instanseg_input)

    // Both masks: Xenium Ranger takes cell boundaries from the cell mask.
    // Nuclei only: Xenium Ranger expands each nucleus by expansion_distance.
    ch_imp_seg_inputs = ch_bundle_path
        .join(INSTANSEG.out.nuclei)
        .join(INSTANSEG.out.cells)
        .combine(ch_coordinate_space)
        .map { meta, bundle, nuclei, cells, coord_space ->
            tuple(
                meta,
                bundle,
                [],
                [],
                nuclei,
                nucleus_segmentation_only ? [] : cells,
                [],
                coord_space,
                expansion_distance,
            )
        }

    XENIUMRANGER_IMPORTSEGMENTATION(ch_imp_seg_inputs)

    emit:
    coordinate_space = ch_coordinate_space                       // channel: [ ["pixels"] ]
    redefined_bundle = XENIUMRANGER_IMPORTSEGMENTATION.out.outs  // channel: [ val(meta), ["redefined-xenium-bundle"] ]
}
