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
  SmithWatermanValidator validator_;
  SmithWatermanValidator& validatorMaster_;
  const std::uint32_t     readLength_;
  WorkQueue&              workQueue_;
  std::mutex&              m_;

  const char* validateOne(const std::vector<std::uint8_t>& block, const char* pCigarEnd);
  void        validateBlock(const std::vector<std::uint8_t>& block);

public:
  ValidatorThread(SmithWatermanValidator& validator, const std::uint32_t readLength, WorkQueue& workQueue, std::mutex& m)
    : validator_(validator), validatorMaster_(validator), readLength_(readLength), workQueue_(workQueue), m_(m)
  {
  }
  void operator()() { runValidator(); }
  void runValidator();
};

}  // namespace simulation
}  // namespace dragenos
