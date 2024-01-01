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
#include "simulation/ReadPackager.hpp"

namespace dragenos {
namespace simulation {

void ReadPackager::flush(){
  workQueue_.acceptBlock(block_, true);
  block_.clear();
  readsIn_ = 0;
}
void ReadPackager::operator()(
      const std::vector<unsigned char>&           query,
      const reference::HashtableConfig::Sequence& contig,
      const std::uint64_t                         refPos,
      const std::uint32_t                         tLen,
      const std::string&                          cigar)
  {
    block_.insert(block_.end(), query.begin(),query.end());
    block_.insert(block_.end(), reinterpret_cast<const unsigned char *>(&contig) , reinterpret_cast<const unsigned char *>(&contig + 1));
    block_.insert(block_.end(), reinterpret_cast<const unsigned char *>(&refPos) , reinterpret_cast<const unsigned char *>(&refPos + 1));
    block_.insert(block_.end(), reinterpret_cast<const unsigned char *>(&tLen) , reinterpret_cast<const unsigned char *>(&tLen + 1));
    std::uint32_t cigarLen = cigar.length();
    block_.insert(block_.end(),reinterpret_cast<const unsigned char *>(&cigarLen) , reinterpret_cast<const unsigned char *>(&cigarLen + 1));
    block_.insert(block_.end(), cigar.begin(), cigar.end());
    ++readsIn_;
    if (readsIn_ == readsPerBlock_){
      assert(!block_.empty());
      workQueue_.acceptBlock(block_, false);
      block_.clear();
      readsIn_ = 0;
    }
  }


}  // namespace simulation
}  // namespace dragenos
