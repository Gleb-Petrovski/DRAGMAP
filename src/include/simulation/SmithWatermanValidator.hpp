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
#pragma once

#include "SmithWatermanRunner.hpp"
#include "simulation/CigarComparer.hpp"
#include "simulation/ReadGenerator.hpp"

namespace dragenos {
namespace simulation {

class SmithWatermanValidator : public ReadGenerator::Processor {
  const SmithWatermanRunner  runner_;
  std::vector<std::uint32_t> histogram_;

public:
  SmithWatermanValidator(const reference::ReferenceDir7& referenceDir) : runner_(referenceDir) {}
  ~SmithWatermanValidator()
  {
    for (std::uint32_t i = 0; i < histogram_.size(); i++) {
      std::cerr << i << ',' << histogram_.at(i) << std::endl;
    }
  }
  virtual void operator()(
      const Query&                                query,
      const reference::HashtableConfig::Sequence& contig,
      const std::uint64_t                         refPos,
      const std::uint32_t                         tLen,
      const std::string&                          cigar) override
  {
    CigarComparer c;

    std::string swCigar;
    histogram_.resize(100 + 1);

    swCigar = runner_.runSW(query, contig, refPos, tLen);

    ++histogram_.at(c.compareCigars(cigar, swCigar) * 100 / c.countMatches(cigar));
    //      if (print < 1 && c.compareCigars(cigar, swCigar) == 0){
    //        std::cerr << cigar << std::endl;
    //        std::cerr << swCigar << std::endl;
    //        std::cerr << convertToString(readSeq, referenceDir) << std::endl;
    //        std::cerr << convertToString(refSeq, referenceDir) <<std::endl;
    //        std::cerr << std::endl;
    //        ++print;
    //      }
  }
};

}  // namespace simulation
}  // namespace dragenos
