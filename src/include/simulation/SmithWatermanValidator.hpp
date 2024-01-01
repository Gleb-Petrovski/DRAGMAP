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

class SmithWatermanValidator {
  const SmithWatermanRunner  runner_;
  const std::uint32_t        flankSizeStart_;
  const std::uint32_t        flankSizeEnd_;
  std::ostream&              output_;
  std::vector<std::uint32_t> histogram_;

public:
  SmithWatermanValidator(
      const reference::ReferenceDir7& referenceDir,
      const uint32_t                  flankSizeStart,
      const uint32_t                  flankSizeEnd,
      std::ostream&                   output)
    : runner_(referenceDir), flankSizeStart_(flankSizeStart), flankSizeEnd_(flankSizeEnd), output_(output)
  {
  }
  void printResults();

  void validate(
        const std::uint8_t*                         queryStart,
        const std::uint8_t*                         queryEnd,
        const reference::HashtableConfig::Sequence& contig,
        const std::uint64_t                         refPos,
        const std::uint32_t                         tLen,
        const std::string&                          cigar);
};

}  // namespace simulation
}  // namespace dragenos
