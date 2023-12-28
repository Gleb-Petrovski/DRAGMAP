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


#include "align/Database.hpp"
#include "align/SmithWaterman.hpp"
#include "simulation/SwSimulator.hpp"

namespace dragenos {
namespace simulation {

std::string SwSimulator::convertToString(const std::vector<unsigned char>& seq, const reference::ReferenceDir7 &referenceDir){
  std::string ret;
  for (const auto& b : seq){
          ret += referenceDir.getReferenceSequence().decodeBase(b);
        }
  return ret;
}

void SwSimulator::runSW(const std::vector<unsigned char>& read, const reference::HashtableConfig::Sequence& s, const reference::ReferenceDir7 &referenceDir
    , std::uint64_t refPos, std::uint32_t readLength,std::string& cigar, std::uint32_t tLen)
{
  const short                             match    = 1;
  const short                             mismatch = -1;
  const dragenos::align::SimilarityScores similarityScores(match, mismatch);
  const int                               gapInit   = 2;
  const int                               gapExtend = 1;
  const std::uint32_t bufferSize = 0;
  align::SmithWaterman sw(similarityScores, gapInit, gapExtend);
  std::string result;
  std::uint32_t startBuffer = std::min<std::uint64_t>(bufferSize, refPos);
  std::uint32_t endBuffer = std::min<std::uint64_t>(bufferSize, s.seqLen - refPos);
  dragenos::align::Database seq;
  referenceDir.getReferenceSequence().getBases(s.seqStart + refPos - startBuffer ,s.seqStart + refPos + readLength + tLen + endBuffer, seq);
  sw.align(&read.front(), &read.front() + read.size(), &seq.front() ,&seq.front() + seq.size(), 10, 2, false, result);
  while (result.at(0) == 'N'){
    refPos++;
    result.erase(result.begin());
  }

//  referenceDir.getReferenceSequence().getBases(s.seqStart + refPos - startBuffer ,s.seqStart + refPos + readLength + tLen + endBuffer, seq);
//

  if (result != cigar)
  {
    std::cerr << convertToString(seq, referenceDir) << std::endl;
        std::cerr << convertToString(read, referenceDir) << std::endl;
        std::cerr << result << std::endl;
        std::cerr << cigar << std::endl;
        std::cerr << std::endl;
  }
}


}  // namespace simulation
}  // namespace dragenos


