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
#include <condition_variable>
#include <mutex>

#include "simulation/SmithWatermanValidator.hpp"

namespace dragenos {
namespace simulation {

class WorkQueue {
  const std::uint32_t maxCount_;
  std::vector<std::uint8_t> block_;
  std::list<std::vector<uint8_t>> queue_;
  std::uint32_t count_ = 0;
  std::mutex                m_;
  std::condition_variable   cv_;
  bool                      ready_     = false;
  bool                      lastBlock_ = false;

public:
  WorkQueue(const std::uint32_t maxCount):maxCount_(maxCount) {}
  void acceptBlock(std::vector<std::uint8_t>& block, bool lastBlock);
  bool getBlock(std::vector<std::uint8_t>& block);
};

}  // namespace simulation
}  // namespace dragenos
