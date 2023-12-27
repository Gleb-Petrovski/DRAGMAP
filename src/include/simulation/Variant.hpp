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
#include <iostream>

namespace dragenos{
namespace simulation{

struct Variant
{
  std::uint64_t refPos_;
  std::uint32_t refLen_;
  std::vector<unsigned char> seq_;
  std::uint64_t refEnd() const
  {
    return refPos_ + refLen_;
  }
  friend std::ostream& operator << (std::ostream& os, const Variant& v)
  {
    return os << "Variant(" << v.refPos_ << "rp " << v.refLen_ << "rl ";
  }

};
typedef std::vector<Variant> Variants;


}
}
