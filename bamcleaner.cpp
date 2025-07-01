####
## cmake   -DSEQ_LIB_INCLUDE_DIR=$HOME/git/svaba/SeqLib   -DSEQ_LIB_LIBRARY=$HOME/git/svaba/build/SeqLib/libseqlib.a   -DHTSLIB_INCLUDE_DIR=/usr/include   -DHTSLIB_LIBRARY=/usr/lib/x86_64-linux-gnu/libhts.so   ..
###


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>
#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"

#include <tsl/robin_set.h>

void print_help(const char* prog) {
  std::cerr << "Usage: " << prog << " <input.bam> <regions.bed> <output.bam>\n"
	    << "  BamCleaner identifies reads in problematic areas\n"
	    << "  and removes read pairs from those regions in a two-pass approach.\n"
	    << "  <regions.bed> is a BED file listing problematic intervals.\n";
}

int main(int argc, char* argv[]) {
  if (argc != 4 || std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help") {
    print_help(argv[0]);
    return (argc == 4) ? 0 : 1;
  }
  
  const std::string inBam   = argv[1];
  const std::string bedFile = argv[2];
  const std::string outBam  = argv[3];
  
  
  // Open BAM reader
  SeqLib::BamReader reader;
  if (!reader.Open(inBam)) {
    std::cerr << "Error: could not open BAM file '" << inBam << "'\n";
    return 1;
  }

  // Load problematic regions from BED
  SeqLib::GRC badRegions;
  std::ifstream bed(bedFile);
  if (!bed) {
    std::cerr << "Error: could not open BED file '" << bedFile << "'\n";
    return 1;
  }
  assert(badRegions.ReadBED(bedFile, reader.Header()));
  badRegions.MergeOverlappingIntervals();
  badRegions.CreateTreeMap();

  std::cerr << "...read " << SeqLib::AddCommas(badRegions.size()) <<
    " BED regions" << std::endl;
  
  // Prepare BAM writer with same header and index mode
  // SeqLib::BamWriter writer;
  // writer.SetHeader(reader.Header());
  // if (!writer.Open(outBam, reader.Header())) {
  //   std::cerr << "Error: could not open output BAM '" << outBam << "'\n";
  //   return 1;
  // }
  
  tsl::robin_set<std::string> qnameset; 
  qnameset.reserve(10'000'000);

  // First pass: flag reads overlapping bad regions
  size_t i = 0;
  size_t added = 0;
  while (auto read = reader.Next()) {
    ++i;
    SeqLib::GenomicRegion rg = read->AsGenomicRegion();
    if (badRegions.CountOverlaps(rg)) {
      const std::string qname = read->Qname();
      if (qnameset.find(qname) != qnameset.end()) {
	++added;
	qnameset.insert(qname);
	if (added % 100'000 == 0)
	  std::cerr
	    << "...added "
	    << std::right               // ensure right-alignment in the field
	    << std::setw(13)            // width = length of "1,000,000,000"
	    << SeqLib::AddCommas(added)
	    << " bad qnames"
	    << std::endl;	  
      }
    }

    if (i % 25'000'000 == 0)
      std::cerr
	<< "...scanned "
	<< std::right               // ensure right-alignment in the field
	<< std::setw(13)            // width = length of "1,000,000,000"
	<< SeqLib::AddCommas(i)
	<< " reads -- " << rg.ChrName(reader.Header())
	<< ":" << std::left << std::setw(12) << SeqLib::AddCommas((rg.pos1 + 1))
	<< std::endl;	  
  }

  // Reset reader for second pass
  //reader.Rewind();
  
  // Second pass: write only reads not in flaggedPairs
  // while (reader.GetNextAlignment(aln)) {
  //   if (flaggedPairs.count(aln.QName()) == 0) {
  //     writer.SaveAlignment(aln);
  //   }
  // }
  
  //writer.Close();
  reader.Close();
  std::cerr << "Finished cleaning BAM.\n";
  return 0;
}
