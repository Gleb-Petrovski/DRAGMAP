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

#include "simulation/Variant.hpp"
namespace dragenos {
namespace simulation {

class VariantGenerator {
public:
  const int maxVarLen_;
  const int varSpacingTarget_;
  VariantGenerator(const int maxVarLen, const int varSpacingTarget)
    : maxVarLen_(maxVarLen), varSpacingTarget_(varSpacingTarget)
  {
  }
  Variants generateVariants(std::uint32_t beginPos, std::uint32_t endPos);

private:
};

}  // namespace simulation
}  // namespace dragenos
