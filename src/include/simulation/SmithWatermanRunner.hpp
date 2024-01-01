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

#include <vector>
#include "reference/ReferenceDir.hpp"

namespace dragenos {
namespace simulation {
typedef std::vector<unsigned char> Query;
class SmithWatermanRunner {
  const reference::ReferenceDir7& referenceDir_;

public:
  SmithWatermanRunner(const reference::ReferenceDir7& referenceDir) : referenceDir_(referenceDir) {}
  std::string runSW(
      const std::uint8_t*                         queryStart,
      const std::uint8_t*                         queryEnd,
      const reference::HashtableConfig::Sequence& s,
      const std::uint32_t                         start,
      const std::uint32_t                         end,
      std::uint64_t                               refPos,
      std::uint32_t                               tLen) const;

private:
};

}  // namespace simulation
}  // namespace dragenos
