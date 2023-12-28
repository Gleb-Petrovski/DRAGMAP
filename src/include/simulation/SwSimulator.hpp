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

#include <vector>
#include "simulation/SwSimulator.hpp"
#include "reference/ReferenceDir.hpp"

namespace dragenos{
namespace simulation{

class SwSimulator{
public:
  SwSimulator()
  {

  }
  void runSW(const std::vector<unsigned char>& read, const reference::HashtableConfig::Sequence& s, const reference::ReferenceDir7 &referenceDir
      , std::uint64_t refPos, std::uint32_t readLength,std::string& cigar, std::uint32_t tLen);
private:
  std::string convertToString(const std::vector<unsigned char>& seq, const reference::ReferenceDir7 &referenceDir);
};

}
}
