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

#include "simulation/WorkQueue.hpp"
#include <assert.h>

namespace dragenos {
namespace simulation {

bool WorkQueue::getBlock(std::vector<std::uint8_t>& block)
{
  bool             ret = false;
  std::unique_lock lk(m_);
  cv_.wait(lk, [this] { return !queue_.empty() || lastBlock_ ; });
  //assert(true == ready_);
  if (!queue_.empty()) {
    queue_.front().swap(block);
    queue_.pop_front();
    --count_;
    ret = true;
  }
  // Manual unlocking is done before notifying, to avoid waking up
  // the waiting thread only to block again (see notify_one for details)
  lk.unlock();
  cv_.notify_all();

  return ret;
}

void WorkQueue::acceptBlock(std::vector<std::uint8_t>& block, bool lastBlock)
{
  std::unique_lock lk(m_);
  cv_.wait(lk, [this] { return queue_.size() <= maxCount_; });
  assert(queue_.size() <= maxCount_);
  queue_.resize(++count_);
  queue_.back().swap(block);
  lastBlock_ = lastBlock;
  // Manual unlocking is done before notifying, to avoid waking up
  // the waiting thread only to block again (see notify_one for details)
  lk.unlock();
  cv_.notify_one();
}

}  // namespace simulation
}  // namespace dragenos
