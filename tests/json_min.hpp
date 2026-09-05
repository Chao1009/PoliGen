// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef LIPOLGEN_TESTS_JSON_MIN_HPP
#define LIPOLGEN_TESTS_JSON_MIN_HPP

// A ~200-line JSON reader for the P1 reference tables.  Deliberately not a
// dependency: the tests must build with nothing but doctest on the include
// path.  Supports the whole of the JSON the dump script emits -- objects,
// arrays, doubles (parsed with strtod, so the shortest-round-trip decimals
// Python writes come back as the same bits), strings with the standard
// escapes, true/false/null.

#include <cctype>
#include <cstdlib>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace jsonmin {

class Value;
using Object = std::map<std::string, Value>;
using Array = std::vector<Value>;

class Value {
 public:
  enum class Type { kNull, kBool, kNumber, kString, kArray, kObject };

  Value() : type_(Type::kNull) {}

  Type type() const { return type_; }
  bool is_null() const { return type_ == Type::kNull; }
  bool is_number() const { return type_ == Type::kNumber; }
  bool is_array() const { return type_ == Type::kArray; }
  bool is_object() const { return type_ == Type::kObject; }

  double num() const {
    if (type_ != Type::kNumber) throw std::runtime_error("json: not a number");
    return number_;
  }
  bool boolean() const {
    if (type_ != Type::kBool) throw std::runtime_error("json: not a bool");
    return bool_;
  }
  const std::string& str() const {
    if (type_ != Type::kString) throw std::runtime_error("json: not a string");
    return string_;
  }
  const Array& arr() const {
    if (type_ != Type::kArray) throw std::runtime_error("json: not an array");
    return array_;
  }
  const Object& obj() const {
    if (type_ != Type::kObject) throw std::runtime_error("json: not an object");
    return object_;
  }

  std::size_t size() const {
    if (type_ == Type::kArray) return array_.size();
    if (type_ == Type::kObject) return object_.size();
    return 0;
  }

  const Value& operator[](std::size_t i) const {
    if (type_ != Type::kArray || i >= array_.size()) {
      throw std::runtime_error("json: array index out of range");
    }
    return array_[i];
  }

  const Value& operator[](const std::string& key) const {
    if (type_ != Type::kObject) throw std::runtime_error("json: not an object");
    const auto it = object_.find(key);
    if (it == object_.end()) throw std::runtime_error("json: no key " + key);
    return it->second;
  }

  bool has(const std::string& key) const {
    return type_ == Type::kObject && object_.count(key) != 0;
  }

  /// Flatten a (possibly nested) array of numbers into a vector.
  std::vector<double> flat() const {
    std::vector<double> out;
    flat_into(out);
    return out;
  }

  friend class Parser;

 private:
  void flat_into(std::vector<double>& out) const {
    if (type_ == Type::kNumber) {
      out.push_back(number_);
    } else if (type_ == Type::kArray) {
      for (const Value& v : array_) v.flat_into(out);
    } else {
      throw std::runtime_error("json: flat() on a non-numeric value");
    }
  }

  Type type_;
  bool bool_ = false;
  double number_ = 0.0;
  std::string string_;
  Array array_;
  Object object_;
};

class Parser {
 public:
  explicit Parser(const std::string& text) : s_(text), i_(0) {}

  Value parse() {
    skip();
    Value v = value();
    skip();
    return v;
  }

 private:
  [[noreturn]] void fail(const std::string& what) const {
    std::ostringstream os;
    os << "json: " << what << " at offset " << i_;
    throw std::runtime_error(os.str());
  }

  void skip() {
    while (i_ < s_.size() && (s_[i_] == ' ' || s_[i_] == '\t' || s_[i_] == '\n'
                              || s_[i_] == '\r')) {
      ++i_;
    }
  }

  bool literal(const char* lit) {
    const std::size_t n = std::char_traits<char>::length(lit);
    if (s_.compare(i_, n, lit) == 0) {
      i_ += n;
      return true;
    }
    return false;
  }

  Value value() {
    if (i_ >= s_.size()) fail("unexpected end of input");
    const char c = s_[i_];
    if (c == '{') return object();
    if (c == '[') return array();
    if (c == '"') {
      Value v;
      v.type_ = Value::Type::kString;
      v.string_ = string();
      return v;
    }
    if (literal("true") || literal("True")) {
      Value v;
      v.type_ = Value::Type::kBool;
      v.bool_ = true;
      return v;
    }
    if (literal("false") || literal("False")) {
      Value v;
      v.type_ = Value::Type::kBool;
      v.bool_ = false;
      return v;
    }
    if (literal("null")) return Value();
    return number();
  }

  Value object() {
    Value v;
    v.type_ = Value::Type::kObject;
    ++i_;  // '{'
    skip();
    if (i_ < s_.size() && s_[i_] == '}') {
      ++i_;
      return v;
    }
    for (;;) {
      skip();
      if (i_ >= s_.size() || s_[i_] != '"') fail("expected a key");
      const std::string key = string();
      skip();
      if (i_ >= s_.size() || s_[i_] != ':') fail("expected ':'");
      ++i_;
      skip();
      v.object_.emplace(key, value());
      skip();
      if (i_ < s_.size() && s_[i_] == ',') {
        ++i_;
        continue;
      }
      if (i_ < s_.size() && s_[i_] == '}') {
        ++i_;
        return v;
      }
      fail("expected ',' or '}'");
    }
  }

  Value array() {
    Value v;
    v.type_ = Value::Type::kArray;
    ++i_;  // '['
    skip();
    if (i_ < s_.size() && s_[i_] == ']') {
      ++i_;
      return v;
    }
    for (;;) {
      skip();
      v.array_.push_back(value());
      skip();
      if (i_ < s_.size() && s_[i_] == ',') {
        ++i_;
        continue;
      }
      if (i_ < s_.size() && s_[i_] == ']') {
        ++i_;
        return v;
      }
      fail("expected ',' or ']'");
    }
  }

  std::string string() {
    ++i_;  // opening quote
    std::string out;
    while (i_ < s_.size() && s_[i_] != '"') {
      char c = s_[i_++];
      if (c != '\\') {
        out.push_back(c);
        continue;
      }
      if (i_ >= s_.size()) fail("unterminated escape");
      const char e = s_[i_++];
      switch (e) {
        case '"': out.push_back('"'); break;
        case '\\': out.push_back('\\'); break;
        case '/': out.push_back('/'); break;
        case 'b': out.push_back('\b'); break;
        case 'f': out.push_back('\f'); break;
        case 'n': out.push_back('\n'); break;
        case 'r': out.push_back('\r'); break;
        case 't': out.push_back('\t'); break;
        case 'u': {
          if (i_ + 4 > s_.size()) fail("short \\u escape");
          const std::string hex = s_.substr(i_, 4);
          i_ += 4;
          const long cp = std::strtol(hex.c_str(), nullptr, 16);
          if (cp < 0x80) {
            out.push_back(static_cast<char>(cp));
          } else if (cp < 0x800) {
            out.push_back(static_cast<char>(0xC0 | (cp >> 6)));
            out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
          } else {
            out.push_back(static_cast<char>(0xE0 | (cp >> 12)));
            out.push_back(static_cast<char>(0x80 | ((cp >> 6) & 0x3F)));
            out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
          }
          break;
        }
        default: fail("unknown escape");
      }
    }
    if (i_ >= s_.size()) fail("unterminated string");
    ++i_;  // closing quote
    return out;
  }

  Value number() {
    const char* start = s_.c_str() + i_;
    char* end = nullptr;
    const double d = std::strtod(start, &end);
    if (end == start) fail("expected a number");
    i_ += static_cast<std::size_t>(end - start);
    Value v;
    v.type_ = Value::Type::kNumber;
    v.number_ = d;
    return v;
  }

  const std::string& s_;
  std::size_t i_;
};

inline Value parse(const std::string& text) { return Parser(text).parse(); }

/// Read and parse a file; returns false (leaving `out` untouched) when the
/// file does not exist, so a test can skip gracefully.
inline bool load_file(const std::string& path, Value& out) {
  std::ifstream in(path.c_str());
  if (!in) return false;
  std::ostringstream buf;
  buf << in.rdbuf();
  const std::string text = buf.str();
  if (text.find_first_not_of(" \t\r\n") == std::string::npos) return false;
  out = parse(text);
  return true;
}

}  // namespace jsonmin

#endif  // LIPOLGEN_TESTS_JSON_MIN_HPP
