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
#include "simulation/WorkQueue.hpp"

namespace dragenos {
namespace simulation {

class ValidatorThread {
  SmithWatermanValidator& validator_;
  const std::uint32_t     readLength_;
  WorkQueue&              workQueue_;

  const char* validateOne(const std::vector<std::uint8_t>& block, const char* pCigarEnd);
  void        validateBlock(const std::vector<std::uint8_t>& block);

public:
  ValidatorThread(SmithWatermanValidator& validator, const std::uint32_t readLength, WorkQueue& workQueue)
    : validator_(validator), readLength_(readLength), workQueue_(workQueue)
  {
  }
  void operator()() { runValidator(); }
  void runValidator();
};

}  // namespace simulation
}  // namespace dragenos
