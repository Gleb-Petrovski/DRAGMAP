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

namespace dragenos {
namespace simulation {

class CigarComparer {
public:
  CigarComparer() {}
  void                 unitTest();
  static std::uint32_t compareCigars(const std::string& cigar1, const std::string& cigar2);
  static std::uint32_t countMatches(const std::string& cigar);

private:
};

}  // namespace simulation
}  // namespace dragenos
