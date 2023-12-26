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
#include <vector>

#include "simulation/VariantGenerator.hpp"

namespace dragenos {
namespace simulation {



Variants VariantGenerator::generateVariants(std::uint32_t beginPos, std::uint32_t endPos)
{
  std::uint32_t cPos = beginPos;
  const std::vector<char> bases = {'A','C','T','G'};
  Variants vars;
  while (cPos < endPos){

    const int randomValue = std::rand();
    Variant v;


    v.refPos_ =  cPos;
    if (randomValue %3 == 0){
      v.refLen_ = 0;
      for (int i = 0; i < randomValue % maxVarLen_ + 1; i++){
        v.seq_ += bases.at(std::rand()%4);
      }
    }
    else if(randomValue %3 == 1)
    {
      v.refLen_ = 1;
      v.seq_ = bases.at(randomValue%4);
    }
    else
    {
      v.refLen_ = randomValue % maxVarLen_ + 1;
      v.seq_ = "";
    }
    vars.push_back(v);
    cPos += v.refLen_ + randomValue % varSpacingTarget_ + 1;
  }
  return vars;
}


}  // namespace simulation
}  // namespace dragenos


