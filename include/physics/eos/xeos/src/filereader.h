/*~-------------------------------------------------------------------~*
 * Copyright (c) 2018 Los Alamos National Security, LLC
 * All rights reserved.
 *~-------------------------------------------------------------------~*/

/**
 * @file filereader.h
 * @author Oleg Korobkin
 * @date 2018-01-06
 * @brief FileReader class for reading various file formats
 */
#ifndef XEOS_FILEREADER_H_
#define XEOS_FILEREADER_H_

#include <string>
#include <fstream>
#include <iostream>

namespace xeos {

class AsciiFileReader {
 private:
  std::ifstream fs;
  std::string   filename;

 public:
  // default constructor: forbid constructor without a filename!
  AsciiFileReader()
      : AsciiFileReader("") {}

  AsciiFileReader(const char *_fn)    { filename.assign(_fn); }

  void SetFilename(const char* _fn)       { filename.assign(_fn); }
  void SetFilename(const std::string _fn) { filename.assign(_fn); }
  std::string GetFilename() { return filename; }

  void Open() {
    if (!fs.is_open())
      fs.open(filename);
    else
      Rewind(); // if file is already open, rewind it
    assert (!fs.fail());
  }

  void Rewind() {
    assert (fs.is_open());
    fs.clear();
    fs.seekg(0);
    assert (!fs.fail());
  }

  void SkipHeader(const int num_lines) {
    std::string str;
    for (int i=0;i<num_lines;i++) {
      assert(!fs.eof());
      std::getline(fs,str);
    }
  }

  /** skips arbitrary number of lines which start with '#'
   *  return the number of lines skipped
   */
  int SkipHashHeader() {
    std::string str;
    std::streampos oldpos;
    int num_lines_header = 0;
    while (std::getline(fs,str)) {
      if (str[0] != '#') {
        fs.seekg(oldpos); // roll back one line
        break;
      }
      oldpos = fs.tellg();
      ++num_lines_header;
    }
    return num_lines_header;
  }

  int NumLines() {
    Open();
    int num_lines = 0;
    std::string str;
    while (std::getline(fs,str))
      ++num_lines;
    return num_lines;
  }

  /**
   * Reade single line in the file, splits it into fields
   * and fills the 'fields' array with values which have
   * been read. Returns total number of fields parsed
   */
  int ReadFields(int num_fields, double* fields) {
    std::string str;
    if (!std::getline(fs, str))
      return 0; // file ended

    std::stringstream ss(str);
    for (int i=0; i<num_fields; ++i) {
      ss >> fields[i];
      if (ss.fail())
        return i; // return less fields than expected
    }
    return num_fields;
  }

  void Close() {
    if (fs.is_open())
      fs.close();
  }
}; // class FileReader

} // namespace
#endif // XEOS_FILEREADER_H_
