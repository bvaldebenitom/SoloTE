import pysam
import sys

inputfile = sys.argv[1]
samfile = pysam.AlignmentFile(inputfile,"rb")

te_bed_file = sys.argv[2]
te_bed_iterator = pysam.tabix_iterator(open(te_bed_file,"r"),pysam.asBed())

outfilename = sys.argv[3]
outfile = pysam.AlignmentFile(outfilename, "wb",template=samfile)

min_overlap = int(sys.argv[4])

for te in te_bed_iterator:
    te_sequence=te[0]
    te_start=int(te[1])
    te_end=int(te[2])
    te_locusname=te[3]
    te_name = (te_locusname.split("|"))[3]
    sam_iterator=samfile.fetch(te_sequence,te_start,te_end)
    for sam_record in sam_iterator:
        if "N" in sam_record.cigarstring:
            length=sam_record.reference_end-sam_record.reference_start+1
            cigar_start = sam_record.reference_start

            for cigartuple in sam_record.cigartuples:
                cigar_op = int(cigartuple[0])
                cigar_op_length = cigartuple[1]
                cigar_end = cigar_start+cigar_op_length

                #if cigar op = M, check overlap with TE
                if cigar_op==0 :
                    if te_end >= cigar_start and te_start <= cigar_end:
                        intersection_start = max(te_start,cigar_start)
                        intersection_end = min(te_end,cigar_end)
                        intersection_length = intersection_end - intersection_start
                        if intersection_length < min_overlap:
                            continue

                        sam_record.set_tag("GX",te_locusname)
                        if sam_record.mapping_quality == 255:
                            sam_record.set_tag("GN","SoloTE|"+te_locusname)
                        else:
                            sam_record.set_tag("GN","SoloTE|"+te_name)
                        outfile.write(sam_record)
                        break

                cigar_start = cigar_end
        else:
            cigar_start = sam_record.reference_start
            cigar_end = sam_record.reference_end
            if te_end >= cigar_start and te_start <= cigar_end:
                 intersection_start = max(te_start,cigar_start)
                 intersection_end = min(te_end,cigar_end)
                 intersection_length = intersection_end - intersection_start
                 if intersection_length < min_overlap:
                      continue
#                 print(str(te_start)+"\t"+str(te_end)+"\t"+str(cigar_start)+"\t"+str(cigar_end)+"\t"+str(intersection_start)+"\t"+str(intersection_end)+"\t"+str(intersection_length)+"\n")
                 sam_record.set_tag("GX",te_locusname)
                 if sam_record.mapping_quality == 255:
                      sam_record.set_tag("GN","SoloTE|"+te_locusname)
                 else:
                      sam_record.set_tag("GN","SoloTE|"+te_name)
                 outfile.write(sam_record)


