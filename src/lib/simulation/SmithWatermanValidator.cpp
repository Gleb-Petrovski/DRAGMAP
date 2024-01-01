/**
 ** DRAGEN Open Source Software
 ** Copyright (c) 2019-2020 Illumina, Inc.
 ** All rights reserved.
 **
 ** This software is provided under the terms and conditions of the
 ** GNU GENERAL PUBLIC LICENSE Version 3
 **
 ** You should have received a copy of the GNU GENERAL PUBLIC LICENSE Version 3
 ** along with this program. If not, see
 ** <https://github.com/illumina/licenses/>.
 **
 **/



#include <assert.h>
#include "simulation/SmithWatermanValidator.hpp"

namespace dragenos {
namespace simulation {


void SmithWatermanValidator::printResults()
  {
    for (std::uint32_t i = 0; i < histogram_.size(); i++) {
      output_ << i << ',' << histogram_.at(i) << '\n';
    }
  }

void SmithWatermanValidator::validate(
      const std::uint8_t*                         queryStart,
      const std::uint8_t*                         queryEnd,
      const reference::HashtableConfig::Sequence& contig,
      const std::uint64_t                         refPos,
      const std::uint32_t                         tLen,
      const std::string&                          cigar)
  {
    CigarComparer c;
    std::string   swCigar;
    histogram_.resize(100 + 1);

    std::uint32_t start = std::max<std::int64_t>(0, refPos - flankSizeStart_);
    std::uint32_t end   = std::min<std::uint64_t>(contig.seqLen, refPos + flankSizeEnd_ + tLen);

    swCigar = runner_.runSW(queryStart, queryEnd, contig, start, end, refPos, tLen);

    ++histogram_.at(c.compareCigars(cigar, swCigar, refPos, start) * 100 / c.countMatches(cigar));
    //      if (print < 1 && c.compareCigars(cigar, swCigar) == 0){
    //        std::cerr << cigar << std::endl;
    //        std::cerr << swCigar << std::endl;
    //        std::cerr << convertToString(readSeq, referenceDir) << std::endl;
    //        std::cerr << convertToString(refSeq, referenceDir) <<std::endl;
    //        std::cerr << std::endl;
    //        ++print;
    //      }
  }


}  // namespace simulation
}  // namespace dragenos
