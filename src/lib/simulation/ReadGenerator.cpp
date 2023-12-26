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
#include "simulation/ReadGenerator.hpp"


namespace dragenos {
namespace simulation {


std::string ReadGenerator::extractRef(
    const reference::ReferenceDir7 &referenceDir,
    const reference::HashtableConfig::Sequence &s, std::uint64_t refPos,
    std::uint32_t matchLen)
{
  std::string ret;
  dragenos::align::Database seq;
  referenceDir.getReferenceSequence().getBases(s.seqStart + refPos,s.seqStart + refPos + matchLen, seq);
  for (const auto b : seq) {
    ret += referenceDir.getReferenceSequence().decodeBase(b);
  }
  return ret;
}

std::string ReadGenerator::generateSeq(std::uint64_t refPos, const reference::HashtableConfig::Sequence& s
    ,std::uint32_t varIdx,const Variants& vars, const reference::ReferenceDir7& referenceDir
    , std::string& cigar)
{
  std::size_t space = readLength_;
  std::string ret;

  while(space)
  {

    while( varIdx < vars.size() && vars.at(varIdx).refPos_ < refPos)
    {
      ++varIdx;

    }
    if (vars.size() == varIdx)
    {
      break;
    }
    const Variant& v = vars.at(varIdx);
    //std::cerr <<"refPos: "<< refPos << "\t" << v << std::endl;
    std::uint32_t matchLen =std::min<std::uint64_t>(v.refPos_ - refPos, space);
    if (matchLen)
    {
      ret += extractRef(referenceDir, s, refPos, matchLen);
      space -= matchLen;
      refPos += matchLen;
      cigar += std::to_string(matchLen) + "M";
    }
    std::uint32_t insLen = std::min(v.seq_.length(), space);
    if (insLen)
    {
      ret += v.seq_.substr(0,insLen);

      space -= insLen;

      if(!v.refLen_)
      {
        cigar += std::to_string(insLen) + "I";
      }
      else if (!v.seq_.empty())
      {
        assert(v.refLen_ == v.seq_.length());
        cigar += std::to_string(insLen) + "X";
      }
    }
    if (v.seq_.empty() && space)
    {
      cigar += std::to_string(v.refLen_) + "D";
    }
    refPos += v.refLen_;
    ++varIdx;

  }
  std::uint32_t matchLen =std::min<std::uint64_t>(s.seqLen - refPos, space);
  if(matchLen)
  {
    ret += extractRef(referenceDir, s, refPos, matchLen);
    cigar += std::to_string(matchLen) + "M";
  }
  return ret;

}

bool insideDeletion(std::uint64_t refPos, const Variant& v)
{
  return (v.seq_.empty() && v.refPos_ <= refPos && v.refEnd() > refPos);

}

void ReadGenerator::generateReads(const reference::HashtableConfig::Sequence& s, const reference::ReferenceDir7& referenceDir
    , const Variants& vars, const std::string& seqName)
{

  std::uint64_t refPos = 0;
  std::uint32_t varIdx = 0;

  while( refPos < s.seqLen)
  {
    std::string cigar;
    std::string readSeq;

    while( varIdx < vars.size() && vars.at(varIdx).refEnd() <= refPos)
    {
      ++varIdx;
    }
    if (vars.size() == varIdx)
    {
      break;
    }
    if (!insideDeletion(refPos, vars.at(varIdx))){
      readSeq = generateSeq(refPos, s, varIdx, vars, referenceDir, cigar);
      output_.printSam(seqName,refPos,cigar, readSeq);
    }
    refPos += readSpacing_;

  }
}


}  // namespace simulation
}  // namespace dragenos


