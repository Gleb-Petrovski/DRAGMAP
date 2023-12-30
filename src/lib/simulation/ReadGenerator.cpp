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

#include "../../include/simulation/SmithWatermanRunner.hpp"
#include "simulation/CigarComparer.hpp"

namespace dragenos {
namespace simulation {


std::vector<unsigned char> ReadGenerator::extractRef(
    const reference::ReferenceDir7 &referenceDir,
    const reference::HashtableConfig::Sequence &s, std::uint64_t refPos,
    std::uint32_t matchLen)
{
  dragenos::align::Database seq;
  referenceDir.getReferenceSequence().getBases(s.seqStart + refPos,s.seqStart + refPos + matchLen, seq);

//  for (const auto& b : seq){
//         std::cerr << referenceDir.getReferenceSequence().decodeBase(b);
//        }
  return seq;

}

std::uint32_t ReadGenerator::generateSeq(std::uint64_t refPos, const reference::HashtableConfig::Sequence& s
    ,std::uint32_t varIdx,const Variants& vars, const reference::ReferenceDir7& referenceDir
    , std::string& cigar, std::vector<unsigned char>& readSeq)
{
  std::size_t space = readLength_;
  std::vector<unsigned char> seq;
  std::uint32_t tLen=0;
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

      seq = extractRef(referenceDir, s, refPos, matchLen);
      readSeq.insert( readSeq.end(), seq.begin(), seq.end() );

      space -= matchLen;
      refPos += matchLen;

      for(std::uint32_t i = 0; i < matchLen; i++){
        cigar += 'M';
      }
    }

    std::uint32_t insLen = std::min(v.seq_.size(), space);

    if (insLen)
    {

      readSeq.insert( readSeq.end(), v.seq_.begin(), v.seq_.begin() + insLen );
      space -= insLen;

      if(!v.refLen_)
      {
        for(std::uint32_t i = 0; i < insLen; i++){
          cigar += 'I';
        }

      }
      else if (!v.seq_.empty())
      {

        assert(v.refLen_ == v.seq_.size());
        for(std::uint32_t i = 0; i < insLen; i++){
          cigar += 'M';
        }

      }
    }
    if (v.seq_.empty() && space)
    {
      for(std::uint32_t i = 0; i < v.refLen_; i++){
        cigar += 'D';
      }
    }
    refPos += v.refLen_;
    tLen += v.refLen_;
    ++varIdx;

  }
  std::uint32_t matchLen =std::min<std::uint64_t>(s.seqLen - refPos, space);
  if(matchLen)
  {
    seq = extractRef(referenceDir, s, refPos, matchLen);
    readSeq.insert( readSeq.end(), seq.begin(), seq.begin() + matchLen);
    for(std::uint32_t i = 0; i < matchLen; i++){
      cigar += 'M';
    }

  }

  return tLen;

}

bool insideDeletion(std::uint64_t refPos, const Variant& v)
{
  return (v.seq_.empty() && v.refPos_ <= refPos && v.refEnd() > refPos);

}
std::string convertToString(const std::vector<unsigned char>& seq, const reference::ReferenceDir7 &referenceDir){
  std::string ret;
  for (const auto& b : seq){
    ret += referenceDir.getReferenceSequence().decodeBase(b);
  }
  return ret;
}
std::string compressCigar(std::string& cigar)
{
  std::string res;
    int num = 0;
    for (std::uint32_t i = 1; i < cigar.size(); i++)
    {
      if (cigar.at(i) == cigar.at(i-1))
      {
        num++;
      }
      else
      {
        res += std::to_string(num+1) + cigar.at(i-1);
        num = 0;
      }
    }
    res += std::to_string(num +1) + cigar.at(cigar.size()-1);
    return res;
}


void ReadGenerator::generateReads(const reference::HashtableConfig::Sequence& s, const reference::ReferenceDir7& referenceDir
    , const Variants& vars, const std::string& seqName, Processor& proc)
{
  std::uint64_t refPos = 0;
  std::uint32_t varIdx = 0;

  while( refPos < s.seqLen)
  {
    Seq readSeq;
    std::string readSeqString;
    while( varIdx < vars.size() && vars.at(varIdx).refEnd() <= refPos)
    {
      ++varIdx;
    }
    if (vars.size() == varIdx)
    {
      break;
    }
    if (!insideDeletion(refPos, vars.at(varIdx))){
      std::string cigar;
      std::uint32_t tLen = generateSeq(refPos, s, varIdx, vars, referenceDir, cigar, readSeq);
      proc(readSeq, s, refPos, readLength_, tLen, cigar);
      output_.printSam(seqName,refPos,cigar, readSeq);
    }
    refPos += readSpacing_;

  }

}


}  // namespace simulation
}  // namespace dragenos

