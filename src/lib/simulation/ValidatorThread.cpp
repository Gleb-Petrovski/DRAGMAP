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

#include "simulation/ValidatorThread.hpp"
#include <assert.h>

namespace dragenos {
namespace simulation {

const char* ValidatorThread::validateOne(const std::vector<std::uint8_t>& block, const char* pCigarEnd)
{
  const std::uint8_t* pQueryStart = reinterpret_cast<const std::uint8_t*>(pCigarEnd);
  const std::uint8_t* pQueryEnd   = pQueryStart + readLength_;

  const reference::HashtableConfig::Sequence* pContig =
      reinterpret_cast<const reference::HashtableConfig::Sequence*>(pQueryEnd);
  const std::uint64_t* pRefPos     = reinterpret_cast<const std::uint64_t*>(pContig + 1);
  const std::uint32_t* pTLen       = reinterpret_cast<const std::uint32_t*>(pRefPos + 1);
  const std::uint32_t* pCigarLen   = reinterpret_cast<const std::uint32_t*>(pTLen + 1);
  const char*          pCigarStart = reinterpret_cast<const char*>(pCigarLen + 1);
  pCigarEnd                        = reinterpret_cast<const char*>(pCigarStart + *pCigarLen);
  std::string cigar(pCigarStart, pCigarEnd);
  validator_.validate(pQueryStart, pQueryEnd, *pContig, *pRefPos, *pTLen, cigar);
  return pCigarEnd;
}

void ValidatorThread::validateBlock(const std::vector<std::uint8_t>& block)
{
  const char* pCigarEnd = reinterpret_cast<const char*>(&block.front());
  do {
    pCigarEnd = validateOne(block, pCigarEnd);
  } while (std::size_t(std::distance(reinterpret_cast<const char*>(&block.front()), pCigarEnd)) <
           block.size());
  assert(size_t(std::distance(reinterpret_cast<const char*>(&block.front()), pCigarEnd)) == block.size());
}

void ValidatorThread::runValidator()
{
  std::vector<std::uint8_t> block;
  while (workQueue_.getBlock(block)) {
    assert(!block.empty());
    validateBlock(block);
    block.clear();
  }
  if (!block.empty()) {
    validateBlock(block);
  }
}

}  // namespace simulation
}  // namespace dragenos
