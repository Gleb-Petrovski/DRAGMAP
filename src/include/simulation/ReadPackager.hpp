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


#include "simulation/ReadGenerator.hpp"
#include "simulation/WorkQueue.hpp"

namespace dragenos {
namespace simulation {

class ReadPackager : public ReadGenerator::Processor {
  const std::uint32_t readsPerBlock_;
  WorkQueue workQueue_;
  std::vector<unsigned char> block_;
  std::uint32_t readsIn_ = 0;

public:
  ReadPackager(std::uint32_t readsPerBlock, WorkQueue& workQueue):readsPerBlock_(readsPerBlock), workQueue_(workQueue)
  {
  }
  void flush();
  virtual void operator()(
      const std::vector<unsigned char>&           query,
      const reference::HashtableConfig::Sequence& contig,
      const std::uint64_t                         refPos,
      const std::uint32_t                         tLen,
      const std::string&                          cigar) override;

};

}  // namespace simulation
}  // namespace dragenos
