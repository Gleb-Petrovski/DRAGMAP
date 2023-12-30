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

namespace dragenos {
namespace simulation {

class SamPrinter : public ReadGenerator::Processor {
  const reference::ReferenceDir7& referenceDir_;
  std::ostream&                   os_;

public:
  SamPrinter(const reference::ReferenceDir7& referenceDir, std::ostream& os)
    : referenceDir_(referenceDir), os_(os)
  {
  }
  virtual void operator()(
      const ReadGenerator::Seq&                   seq,
      const reference::HashtableConfig::Sequence& contig,
      const std::uint64_t                         refPos,
      const std::uint32_t                         tLen,
      const std::string&                          cigar) override
  {
    const auto& seqNames = referenceDir_.getHashtableConfig().getSequenceNames();
    printSam(seqNames.at(contig.id_), refPos, cigar, tLen, seq);
  }

  void printSam(
      const std::string&        refName,
      const std::uint64_t       refPos,
      const std::string&        cigar,
      const std::uint32_t       tLen,
      const ReadGenerator::Seq& seq);

private:
  void printQName(const std::string& seqName, const std::uint64_t refPos, const std::string& cigar);
  static std::string generateQual(const std::vector<unsigned char>& seq);
};

}  // namespace simulation
}  // namespace dragenos
