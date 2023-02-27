#ifndef __UTILS_TABLE_H
#define __UTILS_TABLE_H

#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ranges>

using std::string;
using std::string_view;
using std::cout;
namespace views = std::ranges::views;
using views::iota;


template<typename T>
concept ColumnType = requires(T a, int i) {
    typename T::value_type;
    typename T::iterator;
    { std::begin(a) };
    { std::end(a) };
    { std::size(a) };
    { a[i] };
    { a.at(i) };
};


namespace utils {

template<ColumnType T>
class Table {
    std::vector<string> column_ids{};
    std::unordered_map<string, T> columns{};

    unsigned int precision_ = 8;
    unsigned int column_width_ = 12;

public:
    using column_type = T;
    using item_type = typename T::value_type;

    // Constructors
    Table();
    Table(string_view id, const T & column);
    template<typename ... Args>
    Table(string_view first_id, const T & first_column, const Args & ... other_cols);

    // Add columns
    Table add_columns(string_view id, const T & column);
    template<typename ... Args>
    Table add_columns(string_view first_id, const T & first_column, const Args & ... other_cols);

    // Column access
    T & operator[](const string & id) { return columns[id]; };
    T & at(const string & id) { return columns.at(id); };
    T & at(const string & id) const { return columns.at(id); };

    // Row access (very rough and inefficient)
    auto row(std::size_t n) const {
        std::vector<item_type> row_{};
        for (auto id : column_ids)
            row_.push_back(columns.at(id).at(n));
        return row_;
    }

    // Access useful info
    auto ncols()     const { return columns.size(); }
    auto nrows()     const { return columns.begin()->second.size(); }
    auto size()      const { return ncols() * nrows(); }
    auto width()     const { return column_width_; }
    auto precision() const { return precision_; }
    const std::vector<string> & ids() const { return column_ids; }

    // Setters
    Table set_width(unsigned int n) { column_width_ = n; return *this; }
    Table set_precision(unsigned int n) { precision_ = n; return *this; }

    // Output (plus an overloaded operator<<)
    void print() const;
    void to_csv(string_view filename) const;

    // helper functions for outputting
    void hrule(std::ostream & output) const;
    void header(std::ostream & output) const;
};


template<ColumnType T>
Table<T>
Table<T>::add_columns(string_view id, const T & column){
    if (!columns.empty())
        if (column.size() != nrows())
            throw std::invalid_argument("Wrong size container");
    column_ids.emplace_back(string(id));
    columns.insert(std::make_pair(id, column));
    return *this;
};


template<ColumnType T>
template<typename ... Args>
Table<T>
Table<T>::add_columns(
    string_view first_id,
    const T & first_column,
    const Args & ... other_cols
) {
    add_columns(first_id, first_column);
    if (sizeof...(other_cols) < 2 || sizeof...(other_cols) % 2 != 0)
        throw std::invalid_argument("Invalid number of arguments");
    add_columns(other_cols...);
    return *this;
}


template<ColumnType T>
Table<T>::Table() : column_ids(), columns() {}


template<ColumnType T>
Table<T>::Table(string_view id, const T & column) {
    add_columns(id, column);
}


template<ColumnType T>
template<typename ... Args>
Table<T>::Table(string_view first_id, const T & first_column, const Args & ... other_cols) {
    add_columns(first_id, first_column, other_cols...);
}

template<ColumnType T>
void
Table<T>::hrule(std::ostream & output) const {
    for (auto i : iota(0u, width() * ncols()))
        output << "-";
    output << "\n";
}

template<ColumnType T>
void
Table<T>::header(std::ostream & output) const {
    for (auto id : ids())
        output << std::setw(width()) << id;
    output << "\n";
}


template<ColumnType T>
std::ostream &
operator<<(std::ostream & output, const Table<T> & table) {
    // Set precision
    output << std::setprecision(table.precision());

    // Print header
    table.hrule(output);
    table.header(output);
    table.hrule(output);

    // Print the elements
    for (auto n : iota(0u, table.nrows())) {
        for (auto item : table.row(n))
            output << std::setw(table.width()) << item;
        output << "\n";
    }
    output << std::endl;

    return output;
}


template<ColumnType T>
void
Table<T>::print() const {
    std::cout << *this;
}


template<ColumnType T>
void
Table<T>::to_csv(string_view filename) const {
    std::cout << "   Writing contents onto file '" << filename << "'\n";

    std::fstream file{string(filename), file.out};

    // Set precision
    file << std::setprecision(this->precision());

    auto inner_ids = column_ids | views::take(ncols()-1);
    auto last_id = column_ids.back();

    // Print headers
    for (auto id : inner_ids) file << id << ",";
    file << last_id << "\n";

    // Print columns
    for (auto n : iota(0u, this->nrows())) {
        auto row = this->row(n);
        auto inner_items = row | views::take(ncols() - 1);
        auto last_item = row.back();
        for (auto item : inner_items)
            file << item << ",";
        file << last_item << "\n";
    }
}

}

#endif
