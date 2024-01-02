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

#include "simulation/ReadPackager.hpp"
#include <assert.h>

namespace dragenos {
namespace simulation {

void ReadPackager::flush()
{
  workQueue_.acceptBlock(block_, true);
  block_.clear();
  readsIn_ = 0;
}
void ReadPackager::operator()(
    const std::vector<uint8_t>&           query,
    const reference::HashtableConfig::Sequence& contig,
    const std::uint64_t                         refPos,
    const std::uint32_t                         tLen,
    const std::string&                          cigar)
{
  block_.insert(block_.end(), query.begin(), query.end());
  block_.insert(
      block_.end(), reinterpret_cast<const uint8_t*>(&contig), reinterpret_cast<const uint8_t*>(&contig + 1));
  block_.insert(
      block_.end(), reinterpret_cast<const uint8_t*>(&refPos), reinterpret_cast<const uint8_t*>(&refPos + 1));
  block_.insert(
      block_.end(), reinterpret_cast<const uint8_t*>(&tLen), reinterpret_cast<const uint8_t*>(&tLen + 1));
  std::uint32_t cigarLen = cigar.length();
  block_.insert(
      block_.end(),
      reinterpret_cast<const uint8_t*>(&cigarLen),
      reinterpret_cast<const uint8_t*>(&cigarLen + 1));
  block_.insert(block_.end(), cigar.begin(), cigar.end());
  ++readsIn_;
  if (readsIn_ == readsPerBlock_) {
    assert(!block_.empty());
    workQueue_.acceptBlock(block_, false);
    block_.clear();
    readsIn_ = 0;
  }
}

}  // namespace simulation
}  // namespace dragenos
