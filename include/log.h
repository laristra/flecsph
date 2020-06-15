/*----------------------------------------------------------------------------*
 * Copyright (c) 2016 Los Alamos National Security, LLC
 * All rights reserved.
 *----------------------------------------------------------------------------*/

#pragma once

/*! @file */

#if defined __GNUC__
#include <cxxabi.h>
#include <execinfo.h>
#endif // __GNUC__

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <time.h>

#include <bitset>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>

#include <mpi.h>

#include <mutex>
#include <sstream>

// FIXME: guards?
#if !defined(_MSC_VER)
#include <sys/time.h>
#include <unistd.h>
#endif

#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

//----------------------------------------------------------------------------//
// Default value for tag bits.
//----------------------------------------------------------------------------//

#ifndef LOG_TAG_BITS
#define LOG_TAG_BITS 64
#endif

//----------------------------------------------------------------------------//
// Set the default strip level. All severity levels that are strictly
// less than the strip level will be stripped.
//
// TRACE 0
// INFO  1
// WARN  2
// ERROR 3
// FATAL 4
//----------------------------------------------------------------------------//

#ifndef LOG_STRIP_LEVEL
#define LOG_STRIP_LEVEL 0
#endif

//----------------------------------------------------------------------------//
// Set color output macros depending on whether or not LOG_COLOR_OUTPUT
// is defined.
//----------------------------------------------------------------------------//

#undef COLOR_BLACK
#undef COLOR_DKGRAY
#undef COLOR_RED
#undef COLOR_LTRED
#undef COLOR_GREEN
#undef COLOR_LTGREEN
#undef COLOR_BROWN
#undef COLOR_YELLOW
#undef COLOR_BLUE
#undef COLOR_LTBLUE
#undef COLOR_PURPLE
#undef COLOR_LTPURPLE
#undef COLOR_CYAN
#undef COLOR_LTCYAN
#undef COLOR_LTGRAY
#undef COLOR_WHITE
#undef COLOR_PLAIN

#undef OUTPUT_BLACK
#undef OUTPUT_DKGRAY
#undef OUTPUT_RED
#undef OUTPUT_LTRED
#undef OUTPUT_GREEN
#undef OUTPUT_LTGREEN
#undef OUTPUT_BROWN
#undef OUTPUT_YELLOW
#undef OUTPUT_BLUE
#undef OUTPUT_LTBLUE
#undef OUTPUT_PURPLE
#undef OUTPUT_LTPURPLE
#undef OUTPUT_CYAN
#undef OUTPUT_LTCYAN
#undef OUTPUT_LTGRAY
#undef OUTPUT_WHITE

#define COLOR_BLACK "\033[0;30m"
#define COLOR_DKGRAY "\033[1;30m"
#define COLOR_RED "\033[0;31m"
#define COLOR_LTRED "\033[1;31m"
#define COLOR_GREEN "\033[0;32m"
#define COLOR_LTGREEN "\033[1;32m"
#define COLOR_BROWN "\033[0;33m"
#define COLOR_YELLOW "\033[1;33m"
#define COLOR_BLUE "\033[0;34m"
#define COLOR_LTBLUE "\033[1;34m"
#define COLOR_PURPLE "\033[0;35m"
#define COLOR_LTPURPLE "\033[1;35m"
#define COLOR_CYAN "\033[0;36m"
#define COLOR_LTCYAN "\033[1;36m"
#define COLOR_LTGRAY "\033[0;37m"
#define COLOR_WHITE "\033[1;37m"
#define COLOR_PLAIN "\033[0m"

#define OUTPUT_BLACK(s) COLOR_BLACK << s << COLOR_PLAIN
#define OUTPUT_DKGRAY(s) COLOR_DKGRAY << s << COLOR_PLAIN
#define OUTPUT_RED(s) COLOR_RED << s << COLOR_PLAIN
#define OUTPUT_LTRED(s) COLOR_LTRED << s << COLOR_PLAIN
#define OUTPUT_GREEN(s) COLOR_GREEN << s << COLOR_PLAIN
#define OUTPUT_LTGREEN(s) COLOR_LTGREEN << s << COLOR_PLAIN
#define OUTPUT_BROWN(s) COLOR_BROWN << s << COLOR_PLAIN
#define OUTPUT_YELLOW(s) COLOR_YELLOW << s << COLOR_PLAIN
#define OUTPUT_BLUE(s) COLOR_BLUE << s << COLOR_PLAIN
#define OUTPUT_LTBLUE(s) COLOR_LTBLUE << s << COLOR_PLAIN
#define OUTPUT_PURPLE(s) COLOR_PURPLE << s << COLOR_PLAIN
#define OUTPUT_LTPURPLE(s) COLOR_LTPURPLE << s << COLOR_PLAIN
#define OUTPUT_CYAN(s) COLOR_CYAN << s << COLOR_PLAIN
#define OUTPUT_LTCYAN(s) COLOR_LTCYAN << s << COLOR_PLAIN
#define OUTPUT_LTGRAY(s) COLOR_LTGRAY << s << COLOR_PLAIN
#define OUTPUT_WHITE(s) COLOR_WHITE << s << COLOR_PLAIN

//----------------------------------------------------------------------------//
// Macro utilities.
//----------------------------------------------------------------------------//

#define _log_util_stringify(s) #s
#define _log_stringify(s) _log_util_stringify(s)
#define _log_concat(a, b) a##b

namespace flecsph {

//----------------------------------------------------------------------------//
// Helper functions.
//----------------------------------------------------------------------------//

inline std::string
log_timestamp(bool underscores = false) {
  char stamp[14];
  time_t t = time(0);
  std::string format = underscores ? "%m%d_%H%M%S" : "%m%d %H:%M:%S";
  strftime(stamp, sizeof(stamp), format.c_str(), localtime(&t));
  return std::string(stamp);
} // log_timestamp

template<char C>
std::string
rstrip(const char * file) {
  std::string tmp(file);
  return tmp.substr(tmp.rfind(C) + 1);
} // rstrip

//----------------------------------------------------------------------------//
// Auxiliary types.
//----------------------------------------------------------------------------//

// Options to configure buffered message packet behavior

#ifndef LOG_MAX_MESSAGE_SIZE
#define LOG_MAX_MESSAGE_SIZE 1024
#endif

#ifndef LOG_MAX_PACKET_BUFFER
#define LOG_MAX_PACKET_BUFFER 1024
#endif

#ifndef LOG_PACKET_FLUSH_INTERVAL
#define LOG_PACKET_FLUSH_INTERVAL 100
#endif

//----------------------------------------------------------------------------//
// Packet type.
//----------------------------------------------------------------------------//

struct packet_t {
  static constexpr size_t sec_bytes = sizeof(time_t);
  static constexpr size_t usec_bytes = sizeof(suseconds_t);

  packet_t(const char * msg = nullptr) {
    timeval stamp;
    if(gettimeofday(&stamp, NULL)) {
      std::cerr << "LOG: call to gettimeofday failed!!! " << __FILE__
                << __LINE__ << std::endl;
      std::exit(1);
    } // if

    strncpy(data_, reinterpret_cast<const char *>(&stamp.tv_sec), sec_bytes);
    strncpy(data_ + sec_bytes, reinterpret_cast<const char *>(&stamp.tv_usec),
      usec_bytes);

    std::ostringstream oss;
    oss << msg;

    strcpy(data_ + sec_bytes + usec_bytes, oss.str().c_str());
  } // packet_t

  time_t const & seconds() const {
    return *reinterpret_cast<time_t const *>(data_);
  } // seconds

  suseconds_t const & useconds() const {
    return *reinterpret_cast<suseconds_t const *>(data_ + sec_bytes);
  } // seconds

  const char * message() {
    return data_ + sec_bytes + usec_bytes;
  } // message

  const char * data() const {
    return data_;
  } // data

  size_t bytes() const {
    return sec_bytes + usec_bytes + LOG_MAX_MESSAGE_SIZE;
  } // bytes

  bool operator<(packet_t const & b) {
    return this->seconds() == b.seconds() ? this->useconds() < b.useconds()
                                          : this->seconds() < b.seconds();
  } // operator <

private:
  char data_[sec_bytes + usec_bytes + LOG_MAX_MESSAGE_SIZE];

}; // packet_t

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

// Forward
inline void flush_packets();

struct mpi_state_t {
  static mpi_state_t & instance() {
    static mpi_state_t s;
    return s;
  } // instance

  void init() {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &size_);

    std::thread flusher(flush_packets);
    instance().flusher_thread().swap(flusher);

    initialized_ = true;
  } // init

  bool initialized() {
    return initialized_;
  }

  int rank() {
    return rank_;
  }
  int size() {
    return size_;
  }

  std::thread & flusher_thread() {
    return flusher_thread_;
  }
  std::mutex & packets_mutex() {
    return packets_mutex_;
  }
  std::vector<packet_t> & packets() {
    return packets_;
  }

  bool run_flusher() {
    return run_flusher_;
  }
  void end_flusher() {
    run_flusher_ = false;
  }

private:
  ~mpi_state_t() {
    if(initialized_) {
      end_flusher();
      flusher_thread_.join();
    } // if
  }

  int rank_;
  int size_;
  std::thread flusher_thread_;
  std::mutex packets_mutex_;
  std::vector<packet_t> packets_;
  bool run_flusher_ = true;
  bool initialized_ = false;

}; // mpi_state_t

///
/// Stream buffer type to allow output to multiple targets
/// a la the tee function.
///

//----------------------------------------------------------------------------//
//! The tee_buffer_t type provides a stream buffer that allows output to
//! multiple targets.
//!
//! @ingroup log
//----------------------------------------------------------------------------//

class tee_buffer_t : public std::streambuf
{
public:
  //--------------------------------------------------------------------------//
  //! The buffer_data_t type is used to hold state and the actual low-level
  //! stream buffer pointer.
  //--------------------------------------------------------------------------//

  struct buffer_data_t {
    bool enabled;
    bool colorized;
    std::streambuf * buffer;
  }; // struct buffer_data_t

  //--------------------------------------------------------------------------//
  //! Add a buffer to which output should be written. This also enables
  //! the buffer,i.e., output will be written to it.
  //--------------------------------------------------------------------------//

  void add_buffer(std::string key, std::streambuf * sb, bool colorized) {
    buffers_[key].enabled = true;
    buffers_[key].buffer = sb;
    buffers_[key].colorized = colorized;
  } // add_buffer

  //--------------------------------------------------------------------------//
  //! Enable a buffer so that output is written to it. This is mainly
  //! for buffers that have been disabled and need to be re-enabled.
  //--------------------------------------------------------------------------//

  bool enable_buffer(std::string key) {
    buffers_[key].enabled = true;
    return buffers_[key].enabled;
  } // enable_buffer

  //--------------------------------------------------------------------------//
  //! Disable a buffer so that output is not written to it.
  //--------------------------------------------------------------------------//

  bool disable_buffer(std::string key) {
    buffers_[key].enabled = false;
    return buffers_[key].enabled;
  } // disable_buffer

protected:
  //--------------------------------------------------------------------------//
  //! Override the overflow method. This streambuf has no buffer, so overflow
  //! happens for every character that is written to the string, allowing
  //! us to write to multiple output streams. This method also detects
  //! colorization strings embedded in the character stream and removes
  //! them from output that is going to non-colorized buffers.
  //!
  //! \param c The character to write. This is passed in as an int so that
  //!          non-characters like EOF can be written to the stream.
  //--------------------------------------------------------------------------//

  virtual int overflow(int c) {
    if(c == EOF) {
      return !EOF;
    }
    else {
      // Get the size before we add the current character
      const size_t tbsize = test_buffer_.size();

      // Buffer the output for now...
      test_buffer_.append(1, char(c)); // takes char

      switch(tbsize) {

        case 0:
          if(c == '\033') {
            // This could be a color string, start buffering
            return c;
          }
          else {
            // No match, go ahead and write the character
            return flush_buffer(all_buffers);
          } // if

        case 1:
          if(c == '[') {
            // This still looks like a color string, keep buffering
            return c;
          }
          else {
            // This is some other kind of escape. Write the
            // buffered output to all buffers.
            return flush_buffer(all_buffers);
          } // if

        case 2:
          if(c == '0' || c == '1') {
            // This still looks like a color string, keep buffering
            return c;
          }
          else {
            // This is some other kind of escape. Write the
            // buffered output to all buffers.
            return flush_buffer(all_buffers);
          } // if

        case 3:
          if(c == ';') {
            // This still looks like a color string, keep buffering
            return c;
          }
          else if(c == 'm') {
            // This is a plain color termination. Write the
            // buffered output to the color buffers.
            return flush_buffer(color_buffers);
          }
          else {
            // This is some other kind of escape. Write the
            // buffered output to all buffers.
            return flush_buffer(all_buffers);
          } // if

        case 4:
          if(c == '3') {
            // This still looks like a color string, keep buffering
            return c;
          }
          else {
            // This is some other kind of escape. Write the
            // buffered output to all buffers.
            return flush_buffer(all_buffers);
          } // if

        case 5:
          if(isdigit(c) && (c - '0') < 8) {
            // This still looks like a color string, keep buffering
            return c;
          }
          else {
            // This is some other kind of escape. Write the
            // buffered output to all buffers.
            return flush_buffer(all_buffers);
          } // if

        case 6:
          if(c == 'm') {
            // This is a color string termination. Write the
            // buffered output to the color buffers.
            return flush_buffer(color_buffers);
          }
          else {
            // This is some other kind of escape. Write the
            // buffered output to all buffers.
            return flush_buffer(all_buffers);
          } // if
      } // switch

      return c;
    } // if
  } // overflow

  //--------------------------------------------------------------------------//
  //! Override the sync method so that we sync all of the output buffers.
  //--------------------------------------------------------------------------//

  virtual int sync() {
    int state = 0;

    for(auto b : buffers_) {
      const int s = b.second.buffer->pubsync();
      state = (state != 0) ? state : s;
    } // for

    // Return -1 if one of the buffers had an error
    return (state == 0) ? 0 : -1;
  } // sync

private:
  // Predicate to select all buffers.
  static bool all_buffers(const buffer_data_t & bd) {
    return bd.enabled;
  } // any_buffer

  // Predicate to select color buffers.
  static bool color_buffers(const buffer_data_t & bd) {
    return bd.enabled && bd.colorized;
  } // any_buffer

  // Flush buffered output to buffers that satisfy the predicate function.
  template<typename P>
  int flush_buffer(P && predicate = all_buffers) {
    int eof = !EOF;

    // Put test buffer characters to each buffer
    for(auto b : buffers_) {
      if(predicate(b.second)) {
        for(auto bc : test_buffer_) {
          const int w = b.second.buffer->sputc(bc);
          eof = (eof == EOF) ? eof : w;
        } // for
      } // if
    } // for

    // Clear the test buffer
    test_buffer_.clear();

    // Return EOF if one of the buffers hit the end
    return eof == EOF ? EOF : !EOF;
  } // flush_buffer

  std::unordered_map<std::string, buffer_data_t> buffers_;
  std::string test_buffer_;

}; // class tee_buffer_t

//----------------------------------------------------------------------------//
//! The tee_stream_t type provides a stream class that writes to multiple
//! output buffers.
//----------------------------------------------------------------------------//

struct tee_stream_t : public std::ostream {

  tee_stream_t() : std::ostream(&tee_) {
    tee_.add_buffer("log", std::clog.rdbuf(), true);
  } // tee_stream_t

  tee_stream_t & operator*() {
    return *this;
  } // operator *

  //--------------------------------------------------------------------------//
  //! Add a new buffer to the output.
  //--------------------------------------------------------------------------//

  void add_buffer(std::string key, std::ostream & s, bool colorized = false) {
    tee_.add_buffer(key, s.rdbuf(), colorized);
  } // add_buffer

  //--------------------------------------------------------------------------//
  //! Enable an existing buffer.
  //!
  //! \param key The string identifier of the streambuf.
  //--------------------------------------------------------------------------//

  bool enable_buffer(std::string key) {
    tee_.enable_buffer(key);
    return true;
  } // enable_buffer

  //--------------------------------------------------------------------------//
  //! Disable an existing buffer.
  //!
  //! \param key The string identifier of the streambuf.
  //--------------------------------------------------------------------------//

  bool disable_buffer(std::string key) {
    tee_.disable_buffer(key);
    return false;
  } // disable_buffer

private:
  tee_buffer_t tee_;

}; // struct tee_stream_t

//----------------------------------------------------------------------------//
// Management type.
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
//! The log_t type provides access to logging parameters and configuration.
//!
//! This type provides access to the underlying logging parameters for
//! configuration and information. The flecsph logging functions provide
//! basic logging with an interface that is similar to Google's GLOG
//! and the Boost logging utilities.
//!
//! @note We may want to consider adopting one of these packages
//! in the future.
//!
//! @ingroup log
//----------------------------------------------------------------------------//

class log_t
{
public:
  /// Copy constructor (disabled)
  log_t(const log_t &) = delete;

  /// Assignment operator (disabled)
  log_t & operator=(const log_t &) = delete;

  ///
  /// Meyer's singleton instance.
  ///
  /// \return The singleton instance of this type.
  ///
  static log_t & instance() {
    static log_t c;
    return c;
  } // instance

  ///
  ///
  ///
  void init(std::string active = "none") {
#if defined(LOG_DEBUG)
    std::cerr << COLOR_LTGRAY << "LOG: initializing runtime" << COLOR_PLAIN
              << std::endl;
#endif

    // Because active tags are specified at runtime, it is
    // necessary to maintain a map of the compile-time registered
    // tag names to the id that they get assigned after the log_t
    // initialization (register_tag). This map will be used to populate
    // the tag_bitset_ for fast runtime comparisons of enabled tag groups.

    // Note: For the time being, the map uses actual strings rather than
    // hashes. We should consider creating a const_string_t type for
    // constexpr string creation.

    // Initialize everything to false. This is the default, i.e., "none".
    tag_bitset_.reset();

    // The default group is always active (unscoped). To avoid
    // output for this tag, make sure to scope all LOG output.
    tag_bitset_.set(0);

    if(active == "all") {
      // Turn on all of the bits for "all".
      tag_bitset_.set();
    }
    else if(active != "none") {
      // Turn on the bits for the selected groups.
      std::istringstream is(active);
      std::string tag;
      while(std::getline(is, tag, ',')) {
        if(tag_map_.find(tag) != tag_map_.end()) {
          tag_bitset_.set(tag_map_[tag]);
        }
        else {
          std::cerr << "LOG WARNING: tag " << tag
                    << " has not been registered. Ignoring this group..."
                    << std::endl;
        } // if
      } // while
    } // if

#if defined(LOG_DEBUG)
    std::cerr << COLOR_LTGRAY << "LOG: active tags (" << active << ")"
              << COLOR_PLAIN << std::endl;
#endif

#if defined(LOG_DEBUG)
    std::cerr << COLOR_LTGRAY << "LOG: initializing mpi state" << COLOR_PLAIN
              << std::endl;
#endif

    mpi_state_t::instance().init();

    initialized_ = true;
  } // init

  ///
  /// Return the tag map.
  ///
  const std::unordered_map<std::string, size_t> & tag_map() {
    return tag_map_;
  } // tag_map

  ///
  /// Return the buffered log stream.
  ///
  std::stringstream & buffer_stream() {
    return buffer_stream_;
  } // stream

  ///
  /// Return the log stream.
  ///
  std::ostream & stream() {
    return *stream_;
  } // stream

  ///
  /// Return the log stream predicated on a boolean.
  /// This method interface will allow us to select between
  /// the actual stream and a null stream.
  ///
  std::ostream & severity_stream(bool active = true) {
    return active ? buffer_stream_ : null_stream_;
  } // stream

  ///
  /// Return a null stream to disable output.
  ///
  std::ostream & null_stream() {
    return null_stream_;
  } // null_stream

  ///
  /// Return the tee stream to allow the user to set configuration options.
  /// FIXME: Need a better interface for this...
  ///
  tee_stream_t & config_stream() {
    return *stream_;
  } // stream

  ///
  /// Return the next tag id.
  ///
  size_t register_tag(const char * tag) {
    // If the tag is already registered, just return the previously
    // assigned id. This allows tags to be registered in headers.
    if(tag_map_.find(tag) != tag_map_.end()) {
      return tag_map_[tag];
    } // if

    const size_t id = ++tag_id_;
    assert(id < LOG_TAG_BITS && "Tag bits overflow! Increase LOG_TAG_BITS");
#if defined(LOG_DEBUG)
    std::cerr << COLOR_LTGRAY << "LOG: registering tag " << tag << ": " << id
              << COLOR_PLAIN << std::endl;
#endif
    tag_map_[tag] = id;
    return id;
  } // next_tag

  ///
  /// Return a reference to the active tag (const version).
  ///
  const size_t & active_tag() const {
    return active_tag_;
  } // active_tag

  ///
  /// Return a reference to the active tag (mutable version).
  ///
  size_t & active_tag() {
    return active_tag_;
  } // active_tag

  bool tag_enabled() {
#if defined(LOG_ENABLE_TAGS)

#if defined(LOG_DEBUG)
    auto active_set = tag_bitset_.test(active_tag_) == 1 ? "true" : "false";
    std::cerr << COLOR_LTGRAY << "LOG: tag " << active_tag_ << " is "
              << active_set << COLOR_PLAIN << std::endl;
#endif

    // If the runtime context hasn't been initialized, return true only
    // if the user has enabled externally-scoped messages.
    if(!initialized_) {
#if defined(LOG_ENABLE_EXTERNAL)
      return true;
#else
      return false;
#endif
    } // if

    return tag_bitset_.test(active_tag_);
#else
    return true;
#endif // LOG_ENABLE_TAGS
  } // tag_enabled

  size_t lookup_tag(const char * tag) {
    if(tag_map_.find(tag) == tag_map_.end()) {
      std::cerr << COLOR_YELLOW << "LOG: !!!WARNING " << tag
                << " has not been registered. Ignoring this group..."
                << COLOR_PLAIN << std::endl;
      return 0;
    } // if

    return tag_map_[tag];
  } // lookup_tag

  bool initialized() {
    return initialized_;
  } // initialized

  int rank() {
    return mpi_state_t::instance().rank();
  } // rank

  int size() {
    return mpi_state_t::instance().size();
  } // rank

private:
  ///
  /// Constructor. This method is hidden because we are a singleton.
  ///
  log_t() : null_stream_(0), tag_id_(0), active_tag_(0) {} // log_t

  ~log_t() {
#if defined(LOG_DEBUG)
    std::cerr << COLOR_LTGRAY << "LOG: log_t destructor" << std::endl;
#endif
  }

  bool initialized_ = false;

  tee_stream_t stream_;
  std::stringstream buffer_stream_;
  std::ostream null_stream_;

  size_t tag_id_;
  size_t active_tag_;
  std::bitset<LOG_TAG_BITS> tag_bitset_;
  std::unordered_map<std::string, size_t> tag_map_;

}; // class log_t

void
flush_packets() {
  while(mpi_state_t::instance().run_flusher()) {
    usleep(LOG_PACKET_FLUSH_INTERVAL);

    {
      std::lock_guard<std::mutex> guard(
        mpi_state_t::instance().packets_mutex());

      if(mpi_state_t::instance().packets().size()) {
        std::sort(mpi_state_t::instance().packets().begin(),
          mpi_state_t::instance().packets().end());

        for(auto & p : mpi_state_t::instance().packets()) {
          log_t::instance().stream() << p.message();
        } // for

        mpi_state_t::instance().packets().clear();
      } // if
    } // scope

  } // while
} // flush_packets

//----------------------------------------------------------------------------//
// Tag scope.
//----------------------------------------------------------------------------//

///
/// \class log_tag_scope_t
/// \brief log_tag_scope_t provides an execution scope for which a given
///        tag id is active.
///
/// This type sets the active tag id to the id passed to the constructor,
/// stashing the current active tag. When the instance goes out of scope,
/// the active tag is reset to the stashed value.
///
struct log_tag_scope_t {
  log_tag_scope_t(size_t tag = 0) : stash_(log_t::instance().active_tag()) {
#if defined(LOG_DEBUG)
    std::cerr << COLOR_LTGRAY << "LOG: activating tag " << tag << COLOR_PLAIN
              << std::endl;
#endif

    // Warn users about externally-scoped messages
    if(!log_t::instance().initialized()) {
      std::cerr
        << COLOR_YELLOW << "LOG: !!!WARNING You cannot use "
        << "tag guards for externally scoped messages!!! "
        << "This message will be active if LOG_ENABLE_EXTERNAL is defined!!!"
        << COLOR_PLAIN << std::endl;
    } // if

    log_t::instance().active_tag() = tag;
  } // log_tag_scope_t

  ~log_tag_scope_t() {
    log_t::instance().active_tag() = stash_;
  } // ~log_tag_scope_t

private:
  size_t stash_;

}; // log_tag_scope_t

// Note that none of the tag interface is thread safe. This will have
// to be fixed in the future. One way to do this would be to use TLS
// for the active tag information.
//
// Another feature that would be nice is if the static size_t definition
// failed with helpful information if the user tries to create a tag
// scope for a tag that hasn't been registered.

// Register a tag group with the runtime (log_t). We need the static
// size_t so that tag scopes can be created quickly during execution.
#define log_register_tag(name)                                                 \
  static size_t name##_log_tag_id =                                            \
    flecsph::log_t::instance().register_tag(_log_stringify(name))

// Lookup the tag id
#define log_tag_lookup(name)                                                   \
  flecsph::log_t::instance().lookup_tag(_log_stringify(name))

// Create a new tag scope.
#define log_tag_guard(name)                                                    \
  flecsph::log_tag_scope_t name##_log_tag_scope__(log_tag_lookup(name))

#define log_tag_map() flecsph::log_t::instance().tag_map()

#define send_to_one(message)                                                   \
                                                                               \
  if(mpi_state_t::instance().initialized()) {                                  \
    packet_t pkt(message);                                                     \
                                                                               \
    packet_t * pkts = mpi_state_t::instance().rank() == 0                      \
                        ? new packet_t[mpi_state_t::instance().size()]         \
                        : nullptr;                                             \
                                                                               \
    MPI_Gather(pkt.data(), pkt.bytes(), MPI_BYTE, pkts, pkt.bytes(), MPI_BYTE, \
      0, MPI_COMM_WORLD);                                                      \
                                                                               \
    if(mpi_state_t::instance().rank() == 0) {                                  \
                                                                               \
      std::lock_guard<std::mutex> guard(                                       \
        mpi_state_t::instance().packets_mutex());                              \
                                                                               \
      for(size_t i{0}; i < mpi_state_t::instance().size(); ++i) {              \
        mpi_state_t::instance().packets().push_back(pkts[i]);                  \
      } /* for */                                                              \
                                                                               \
      delete[] pkts;                                                           \
                                                                               \
    } /* if */                                                                 \
  } /* if */

//----------------------------------------------------------------------------//
// Base type for log messages.
//----------------------------------------------------------------------------//

///
/// Function always returning true. Used for defaults.
///
inline bool
true_state() {
  return true;
} // output_bool

/*!
  The log_message_t type provides a base class for implementing
  formatted logging utilities.
 */
template<typename P>
struct log_message_t {
  /*!
    Constructor.

    @tparam P Predicate function type.

    @param file            The current file (where the log message was
                           created).  In general, this will always use the
                           __FILE__ parameter from the calling macro.
    @param line            The current line (where the log message was called).
                           In general, this will always use the __LINE__
                           parameter from the calling macro.
    @param predicate       The predicate function to determine whether or not
                           the calling runtime should produce output.
    @param can_send_to_one A boolean indicating whether the calling scope
                           can route messages through one rank.
   */
  log_message_t(const char * file,
    int line,
    P && predicate,
    bool can_send_to_one = true)
    : file_(file), line_(line), predicate_(predicate),
      can_send_to_one_(can_send_to_one), clean_color_(false) {
#if defined(LOG_DEBUG)
    std::cerr << COLOR_LTGRAY << "LOG: log_message_t constructor " << file
              << " " << line << COLOR_PLAIN << std::endl;
#endif
  } // log_message_t

  virtual ~log_message_t() {
#if defined(LOG_DEBUG)
    std::cerr << COLOR_LTGRAY << "LOG: log_message_t destructor " << COLOR_PLAIN
              << std::endl;
#endif

    if(can_send_to_one_) {
      send_to_one(log_t::instance().buffer_stream().str().c_str());
    }
    else {
      log_t::instance().stream() << log_t::instance().buffer_stream().str();
    } // if

    log_t::instance().buffer_stream().str(std::string{});
  } // ~log_message_t

  ///
  /// Return the output stream. Override this method to add additional
  // formatting to a particular severity output.
  ///
  virtual std::ostream & stream() {
    return log_t::instance().severity_stream(predicate_());
  } // stream

protected:
  const char * file_;
  int line_;
  P & predicate_;
  bool can_send_to_one_;
  bool clean_color_;

}; // struct log_message_t

//----------------------------------------------------------------------------//
// Convenience macro to define severity levels.
//
// Log message types defined using this macro always use the default
// predicate function, true_state().
//----------------------------------------------------------------------------//

#define log_severity_message_t(severity, P, format)                            \
  struct severity##_log_message_t : public log_message_t<P> {                  \
    severity##_log_message_t(const char * file,                                \
      int line,                                                                \
      P && predicate = true_state,                                             \
      bool can_send_to_one = true)                                             \
      : log_message_t<P>(file, line, predicate, can_send_to_one) {}            \
                                                                               \
    ~severity##_log_message_t() {                                              \
      /* Clean colors from the stream */                                       \
      if(clean_color_) {                                                       \
        log_t::instance().buffer_stream() << COLOR_PLAIN;                      \
      }                                                                        \
    }                                                                          \
                                                                               \
    std::ostream &                                                             \
    stream() override /* This is replaced by the scoped logic */               \
      format                                                                   \
  }

//----------------------------------------------------------------------------//
// Define the insertion style severity levels.
//----------------------------------------------------------------------------//

#define log_message_stamp                                                      \
  log_timestamp() << " " << flecsph::rstrip<'/'>(file_) << ":" << line_

#define log_mpi_stamp " r" << mpi_state_t::instance().rank()

// Trace
log_severity_message_t(trace, decltype(flecsph::true_state), {
  std::ostream & stream = log_t::instance().severity_stream(
    LOG_STRIP_LEVEL < 1 && predicate_() && log_t::instance().tag_enabled());

  {
    stream << OUTPUT_CYAN("[T") << OUTPUT_LTGRAY(log_message_stamp);
    stream << OUTPUT_DKGRAY(log_mpi_stamp);
    stream << OUTPUT_CYAN("] ");
  } // scope

  return stream;
});

// Info
log_severity_message_t(info, decltype(flecsph::true_state), {
  std::ostream & stream = log_t::instance().severity_stream(
    LOG_STRIP_LEVEL < 2 && predicate_() && log_t::instance().tag_enabled());

  {
    stream << OUTPUT_GREEN("[I") << OUTPUT_LTGRAY(log_message_stamp);
    stream << OUTPUT_DKGRAY(log_mpi_stamp);
    stream << OUTPUT_GREEN("] ");
  } // scope

  return stream;
});

// Warn
log_severity_message_t(warn, decltype(flecsph::true_state), {
  std::ostream & stream = log_t::instance().severity_stream(
    LOG_STRIP_LEVEL < 3 && predicate_() && log_t::instance().tag_enabled());

  {
    stream << OUTPUT_BROWN("[W") << OUTPUT_LTGRAY(log_message_stamp);
    stream << OUTPUT_DKGRAY(log_mpi_stamp);
    stream << OUTPUT_BROWN("] ") << COLOR_YELLOW;
  } // scope

  clean_color_ = true;
  return stream;
});

// Error
log_severity_message_t(error, decltype(flecsph::true_state), {
  std::ostream & stream = std::cerr;

  {
    stream << OUTPUT_RED("[E") << OUTPUT_LTGRAY(log_message_stamp);
    stream << OUTPUT_DKGRAY(log_mpi_stamp);
    stream << OUTPUT_RED("] ") << COLOR_LTRED;
  } // scope

  clean_color_ = true;
  return stream;
});

} // namespace flecsph

//----------------------------------------------------------------------------//
// Private macro interface
//----------------------------------------------------------------------------//

//
// Indirection to expand counter name.
//

#define log_counter_varname(str, line) _log_concat(str, line)

//
// Indirection to expand counter name.
//

#define log_counter(str) log_counter_varname(str, __LINE__)

//
// Define a counter name.
//

#define log_counter_name log_counter(counter)

//----------------------------------------------------------------------------//
// Macro Interface
//----------------------------------------------------------------------------//

/*!
  @def log_init(active)

  This call initializes the log runtime with the list of tags specified
  in \em active.

  @param active A const char * or std::string containing the list of
                active tags. Tags should be comma delimited.

  \b Usage
  \code
  int main(int argc, char ** argv) {

     // Fill a string object with the active tags.
     std::string tags{"init,advance,analysis"};

     // Initialize the log runtime with active tags, 'init', 'advance',
     // and 'analysis'.
     log_init(tags);

     return 0;
  } // main
  \endcode

  @ingroup log
 */

#define log_init(active)                                                       \
  /* MACRO IMPLEMENTATION */                                                   \
                                                                               \
  flecsph::log_t::instance().init(active)

/*!
  @def log(severity)

  This handles all of the different logging modes for the insertion
  style logging interface.

  @param severity The severity level of the log entry.

  @note The form "true && ..." is necessary for tertiary argument
        evaluation so that the std::ostream & returned by the stream()
        function can be implicitly converted to an int.

  @b Usage
  @code
  int value{20};

  // Print the value at info severity level
  log(info) << "Value: " << value << std::endl;

  // Print the value at warn severity level
  log(warn) << "Value: " << value << std::endl;
  @endcode

  @ingroup log
 */

#define logm(severity)                                                         \
  /* MACRO IMPLEMENTATION */                                                   \
                                                                               \
  true && flecsph::severity##_log_message_t(__FILE__, __LINE__).stream()

/*!
  @def log_trace(message)

  Method style interface for trace level severity log entries.

  @param message The stream message to be printed.

  @b Usage
  @code
  int value{20};

  // Print the value at trace severity level
  log_trace("Value: " << value);
  @endcode

  @ingroup log
 */

#define log_trace(message)                                                     \
  /* MACRO IMPLEMENTATION */                                                   \
                                                                               \
  flecsph::trace_log_message_t(__FILE__, __LINE__).stream() << message

/*!
  @def log_info(message)

  Method style interface for info level severity log entries.

  @param message The stream message to be printed.

  @b Usage
  @code
  int value{20};

  // Print the value at info severity level
  log_info("Value: " << value);
  @endcode

  @ingroup log
 */

#define log_info(message)                                                      \
  /* MACRO IMPLEMENTATION */                                                   \
                                                                               \
  flecsph::info_log_message_t(__FILE__, __LINE__).stream() << message

/*!
  @def log_warn(message)

  Method style interface for warn level severity log entries.

  @param message The stream message to be printed.

  @b Usage
  @code
  int value{20};

  // Print the value at warn severity level
  log_warn("Value: " << value);
  @endcode

  @ingroup log
 */

#define log_warn(message)                                                      \
  /* MACRO IMPLEMENTATION */                                                   \
                                                                               \
  flecsph::warn_log_message_t(__FILE__, __LINE__).stream() << message

/*!
  @def log_error(message)

  Method style interface for error level severity log entries.

  @param message The stream message to be printed.

  @b Usage
  @code
  int value{20};

  // Print the value at error severity level
  log_error("Value: " << value);
  @endcode

  @ingroup log
 */

#define log_error(message)                                                     \
  /* MACRO IMPLEMENTATION */                                                   \
                                                                               \
  flecsph::error_log_message_t(__FILE__, __LINE__).stream() << message

/*!
  @def log_fatal(message)

  Throw a runtime exception with the provided message.

  @param message The stream message to be printed.

  @note Fatal level severity log entires are not disabled by tags or
        by the ENABLE_LOG or LOG_STRIP_LEVEL build options, i.e.,
        they are always active.

  @b Usage
  @code
  int value{20};

  // Print the value and exit
  log_fatal("Value: " << value);
  @endcode

  @ingroup log
 */

#define log_fatal(message)                                                     \
  /* MACRO IMPLEMENTATION */                                                   \
                                                                               \
  {                                                                            \
    std::stringstream _sstream;                                                \
    _sstream << OUTPUT_LTRED("FATAL ERROR ")                                   \
             << OUTPUT_YELLOW(                                                 \
                  flecsph::rstrip<'/'>(__FILE__) << ":" << __LINE__ << " ")    \
             << OUTPUT_LTRED(message) << std::endl;                            \
    throw std::runtime_error(_sstream.str());                                  \
  } /* scope */

/*!
  @def log_assert(test, message)

  log assertion interface. Assertions allow the developer to catch
  invalid program state. This call will invoke log_fatal if the test
  condition is false.

  @param test    The test condition.
  @param message The stream message to be printed.

  @note Failed assertions are not disabled by tags or
        by the ENABLE_LOG or LOG_STRIP_LEVEL build options, i.e.,
        they are always active.

  @b Usage
  @code
  int value{20};

  // Print the value and exit
  log_assert(value == 20, "invalid value");
  @endcode

  @ingroup log
 */

#define log_assert(test, message)                                              \
  /* MACRO IMPLEMENTATION */                                                   \
                                                                               \
  if(!(test)) {                                                                \
    log_fatal(message);                                                        \
  }
//  !(test) && log_fatal(message)

/*!
  @def log_add_buffer(name, ostream, colorized)

  Add a named stream buffer to the log runtime. Added buffers are enabled
  by default, and can be disabled by calling \ref log_disable_buffer.

  @param name      The name of the output buffer.
  @param ostream   The output stream of type std::ostream.
  @param colorized A boolean indicating whether or not the output to
                   this stream should be colorized.

  @ingroup log
 */

#define log_add_buffer(name, ostream, colorized)                               \
  /* MACRO IMPLEMENTATION */                                                   \
                                                                               \
  flecsph::log_t::instance().config_stream().add_buffer(                       \
    name, ostream, colorized)

/*!
  @def log_enable_buffer(name)

  Enable an output buffer.

  @param name The name of the output stream that was used to add the buffer.

  @ingroup log
 */

#define log_enable_buffer(name)                                                \
  /* MACRO IMPLEMENTATION */                                                   \
                                                                               \
  flecsph::log_t::instance().config_stream().enable_buffer(name)

/*!
  @def log_disable_buffer(name)

  Disable an output buffer.

  @param name The name of the output stream that was used to add the buffer.

  @ingroup log
 */

#define log_disable_buffer(name)                                               \
  /* MACRO IMPLEMENTATION */                                                   \
                                                                               \
  flecsph::log_t::instance().config_stream().disable_buffer(name)

namespace flecsph_log {

/*!
  Enum type to specify output delimiters for containers.

  @ingroup flecsph_log
 */

enum log_delimiters_t : size_t {
  newline,
  space,
  colon,
  semicolon,
  comma
}; // enum log_delimiters_t

} // namespace flecsph_log

// \TODO actually fix warning
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wtautological-compare"
#endif

/*!
  @def log_container(severity, banner, container, delimiter)

  Output the contents of a standard container type. Valid container types
  must implement a forward iterator.

  @param severity  The severity level at which to output the message.
  @param banner    A top-level label for the container output.
  @param container The container to output.
  @param delimiter The output character to use to delimit container
                   entries, e.g., newline, comma, space, etc. Valid
                   delimiters are defined in log_delimiters_t.

  @ingroup flecsph_log
 */

#define log_container(severity, banner, container, delimiter)                  \
  /* MACRO IMPLEMENTATION */                                                   \
                                                                               \
  {                                                                            \
    std::stringstream ss;                                                      \
    char delim =                                                               \
      (delimiter == flecsph_log::newline)                                      \
        ? '\n'                                                                 \
        : (delimiter == flecsph_log::space)                                    \
            ? ' '                                                              \
            : (delimiter == flecsph_log::colon)                                \
                ? ':'                                                          \
                : (delimiter == flecsph_log::semicolon) ? ';' : ',';           \
    ss << banner << (delimiter == flecsph_log::newline ? '\n' : ' ');          \
    size_t entry(0);                                                           \
    for(auto c = container.begin(); c != container.end(); ++c) {               \
      (delimiter == flecsph_log::newline) &&                                   \
        ss << OUTPUT_CYAN("[C") << OUTPUT_LTGRAY(" entry ") << entry++         \
           << OUTPUT_CYAN("]") << std::endl                                    \
           << *c;                                                              \
      (c != --container.end()) && ss << delim;                                 \
    }                                                                          \
    logm(severity) << ss.str() << std::endl;                                   \
  }

#if defined(__clang__)
#pragma clang diagnostic pop
#endif

// Enable MPI
namespace flecsph {

/*!
  The mpi_config_t type provides an interface to MPI runtime state
  information.

  @ingroup flecsph_log
 */

struct mpi_config_t {

  /*!
    Meyer's singleton instance.

    @return The single instance of this type.
   */

  static mpi_config_t & instance() {
    static mpi_config_t m;
    return m;
  } // instance

  /*!
    Return the active rank as a constant reference.
   */

  const size_t & active_rank() const {
    return active_rank_;
  } // active_rank

  /*!
    Return the active rank as a mutable reference.
   */

  size_t & active_rank() {
    return active_rank_;
  } // active_rank

private:
  mpi_config_t() : active_rank_(0) {}

  size_t active_rank_;

}; // struct mpi_config_t

/*!
  Return a boolean indicating whether the current runtime rank matches a
  statically defined value.

  @tparam RANK The static rank to use in the comparison.

  @ingroup flecsph_log
 */

template<size_t RANK>
inline bool
is_static_rank() {
  int part;
  MPI_Comm_rank(MPI_COMM_WORLD, &part);
  return part == RANK;
} // is_static_rank

/*!
  Return a boolean that indicates whether the current runtime rank is active.

  @ingroup flecsph_log
 */

inline bool
is_active_rank() {
  int part;
  MPI_Comm_rank(MPI_COMM_WORLD, &part);
  return part == mpi_config_t::instance().active_rank();
} // is_active_rank

} // namespace flecsph

/*!
  @def log_rank(severity, rank)

  This handles all of the different logging modes for the insertion
  style logging interface.

  @param severity The severity level of the log entry.
  @param rank     The rank for which to output the message stream.

  @note The form "true && ..." is necessary for tertiary argument
        evaluation so that the std::ostream & returned by the stream()
        function can be implicitly converted to an int.

  @b Usage
  @code
  int value{20};

  // Print the value at info severity level on rank 0
  log_rank(info, 0) << "Value: " << value << std::endl;

  // Print the value at warn severity level on rank 1
  log_rank(warn, 1) << "Value: " << value << std::endl;
  @endcode

  @ingroup flecsph_log
 */

#define log_rank(severity, rank)                                               \
  /* MACRO IMPLEMENTATION */                                                   \
                                                                               \
  true && flecsph::severity##_log_message_t(                                   \
            __FILE__, __LINE__, flecsph::is_static_rank<rank>, false)          \
            .stream()

/*!
  @def log_set_output_rank(rank)

  Set the output rank for calls to log_one.

  @param rank The rank for which output will be generated.

  @ingroup flecsph_log
 */

#define log_set_output_rank(rank)                                              \
  /* MACRO IMPLEMENTATION */                                                   \
                                                                               \
  flecsph::mpi_config_t::instance().active_rank() = rank

/*!
  @def log_one(severity)

  This handles all of the different logging modes for the insertion
  style logging interface. This will only output on the rank specified
  by log_set_output_rank.

  @param severity The severity level of the log entry.

  @ingroup flecsph_log
 */

#define log_one(severity)                                                      \
  /* MACRO IMPLEMENTATION */                                                   \
                                                                               \
  true && flecsph::severity##_log_message_t(                                   \
            __FILE__, __LINE__, flecsph::is_active_rank, false)                \
            .stream()

/*!
  @def log_container_rank(severity, banner, container, delimiter, rank)

  Output the contents of a standard container type on the specified
  rank. Valid container types must implement a forward iterator.

  @param severity  The severity level at which to output the message.
  @param banner    A top-level label for the container output.
  @param container The container to output.
  @param delimiter The output character to use to delimit container
                   entries, e.g., newline, comma, space, etc. Valid
                   delimiters are defined in log_delimiters_t.
  @param rank      The rank for which to output the message stream.

  @ingroup flecsph_log
 */

#define log_container_rank(severity, banner, container, delimiter, rank)       \
  /* MACRO IMPLEMENTATION */                                                   \
                                                                               \
  {                                                                            \
    std::stringstream ss;                                                      \
    char delim =                                                               \
      (delimiter == flecsph_log::newline)                                      \
        ? '\n'                                                                 \
        : (delimiter == flecsph_log::space)                                    \
            ? ' '                                                              \
            : (delimiter == flecsph_log::colon)                                \
                ? ':'                                                          \
                : (delimiter == flecsph_log::semicolon) ? ';' : ',';           \
    ss << banner << (delimiter == flecsph_log::newline ? '\n' : ' ');          \
    size_t entry(0);                                                           \
    for(auto c = container.begin(); c != container.end(); ++c) {               \
      (delimiter == flecsph_log::newline) &&                                   \
        ss << OUTPUT_CYAN("[C") << OUTPUT_LTGRAY(" entry ") << entry++         \
           << OUTPUT_CYAN("]") << std::endl;                                   \
      ss << *c;                                                                \
      (c != --container.end()) && ss << delim;                                 \
    }                                                                          \
    log_rank(severity, rank) << ss.str() << std::endl;                         \
  } /* scope */

/*!
  @def log_container_one(severity, banner, container, delimiter)

  Output the contents of a standard container type on the rank
  specified by log_set_output_rank. Valid container types must
  implement a forward iterator.

  @param severity  The severity level at which to output the message.
  @param banner    A top-level label for the container output.
  @param container The container to output.
  @param delimiter The output character to use to delimit container
                   entries, e.g., newline, comma, space, etc. Valid
                   delimiters are defined in log_delimiters_t.

  @ingroup flecsph_log
 */

#define log_container_one(severity, banner, container, delimiter)              \
  /* MACRO IMPLEMENTATION */                                                   \
                                                                               \
  {                                                                            \
    std::stringstream ss;                                                      \
    char delim =                                                               \
      (delimiter == flecsph_log::newline)                                      \
        ? '\n'                                                                 \
        : (delimiter == flecsph_log::space)                                    \
            ? ' '                                                              \
            : (delimiter == flecsph_log::colon)                                \
                ? ':'                                                          \
                : (delimiter == flecsph_log::semicolon) ? ';' : ',';           \
    ss << banner << (delimiter == flecsph_log::newline ? '\n' : ' ');          \
    size_t entry(0);                                                           \
    for(auto c = container.begin(); c != container.end(); ++c) {               \
      (delimiter == flecsph_log::newline) &&                                   \
        ss << OUTPUT_CYAN("[C") << OUTPUT_LTGRAY(" entry ") << entry++         \
           << OUTPUT_CYAN("]") << std::endl;                                   \
      ss << *c;                                                                \
      (c != --container.end()) && ss << delim;                                 \
    }                                                                          \
    log_one(severity) << ss.str() << std::endl;                                \
  } /* scope */
