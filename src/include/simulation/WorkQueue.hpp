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


#include "simulation/SmithWatermanValidator.hpp"

namespace dragenos {
namespace simulation {

class WorkQueue{
  SmithWatermanValidator& validator_;
  std::uint32_t readLength_;
public:
  WorkQueue(SmithWatermanValidator& validator, const std::uint32_t readLength):validator_(validator), readLength_(readLength)
  {
  }
  void acceptBlock(const std::vector<std::uint8_t> block);

};

}  // namespace simulation
}  // namespace dragenos
