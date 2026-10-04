// subworkflows/local/prepare_crams/main.nf
//
// Make every input CRAM readable by GATK 4.5: pass CRAM 3.0 files through and
// rewrite any other version (in practice 3.1) with CRAM_TO_V30. The version is
// read from the 6-byte file definition at the start of each CRAM ("CRAM",
// major, minor), so the check costs nothing at launch.
include { CRAM_TO_V30 } from '../../../modules/local/cram_to_v30/main.nf'

def cramVersion(cram) {
    def b = new byte[6]
    cram.withInputStream { it.read(b) }
    if (new String(b, 0, 4, 'US-ASCII') != 'CRAM') {
        error "${cram} does not start with the CRAM magic bytes; is it a CRAM file?"
    }
    return "${b[4]}.${b[5]}".toString()
}

workflow PREPARE_CRAMS {
    take:
    cram_ch        // channel: [strain, cram, crai]
    reference
    reference_fai

    main:
    branched = cram_ch.branch { strain, cram, crai ->
        v30:   cramVersion(cram) == '3.0'
        other: true
    }
    branched.other.subscribe { strain, cram, crai ->
        log.warn "${strain}: ${cram.name} is CRAM ${cramVersion(cram)}; GATK 4.5 reads only CRAM 3.0, rewriting it with CRAM_TO_V30"
    }
    CRAM_TO_V30(branched.other, reference, reference_fai)

    emit:
    cram = branched.v30.mix(CRAM_TO_V30.out.cram)
}
