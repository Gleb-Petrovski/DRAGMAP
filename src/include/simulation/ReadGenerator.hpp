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
#include "reference/ReferenceDir.hpp"
#include "simulation/Variant.hpp"
#include "simulation/SamPrinter.hpp"

namespace dragenos{
namespace simulation{

class ReadGenerator{
public:
  SamPrinter& output_;
  const int readLength_;
  const int readSpacing_;
  ReadGenerator(SamPrinter& output, const int readLength, const int readSpacing):output_(output)
  , readLength_(readLength), readSpacing_(readSpacing)
  {
  }
  void generateReads(const reference::HashtableConfig::Sequence& s, const reference::ReferenceDir7& referenceDir
      , const Variants& vars,const std::string& seqName);
private:
  std::vector<unsigned char> extractRef(
      const reference::ReferenceDir7 &referenceDir,
      const reference::HashtableConfig::Sequence &s, std::uint64_t refPos,
      std::uint32_t matchLen);
  std::uint32_t generateSeq(std::uint64_t refPos, const reference::HashtableConfig::Sequence& s
      ,std::uint32_t varIdx,const Variants& vars, const reference::ReferenceDir7& referenceDir
      , std::string& cigar, std::vector<unsigned char>& readSeq);
};


}
}
