#ifndef __UTILS_TABLE_H
#define __UTILS_TABLE_H

#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <stdexcept>
#include <ostream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ranges>

using std::string;
using std::string_view;
using std::cout;
using std::setw;
namespace views = std::ranges::views;
using views::iota;

template<typename T>
concept ColumnType = requires(T a, int i) {
    { std::begin(a) };
    { std::end(a) };
    { std::size(a) };
    { a[i] };
};

namespace utils {

template<ColumnType T>
class Table {
    std::vector<string> column_ids{};
    std::unordered_map<string, T> columns{};

    unsigned int precision = 8;
    unsigned int col_width = 12;

public:
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

    const std::vector<string> & head() const { return column_ids; }

    auto get_width() const { return col_width; }
    auto get_precision() const { return precision; }

    Table set_width(unsigned int n) {
        col_width = n;
        return *this;
    }

    Table set_precision(unsigned int n) {
        precision = n;
        return *this;
    }

    auto ncols() const { return columns.size(); }
    auto nrows() const { return columns.begin()->second.size(); }
    auto size()  const { return ncols() * nrows(); }

    // Output
    void print() const;
    void to_csv(string_view filename) const;

    T & operator[](const string & id) {
        return columns[id];
    };

    auto row(std::size_t n) const {
        std::vector<item_type> row_{};
        for (auto id : column_ids)
            row_.push_back(columns.at(id).at(n));
        return row_;
    }

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
Table<T>::add_columns(string_view first_id, const T & first_column, const Args & ... other_cols)
{
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
std::ostream &
operator<<(std::ostream & output, const Table<T> & table) {
    // Set precision
    output << std::setprecision(table.get_precision());

    // Print the headers
    for (auto id : table.head())
        output << setw(table.get_width()) << id;
    output << "\n";

    // Print hrule
    for (auto i : iota(0u, table.get_width() * table.ncols()))
        output << "-";
    output << "\n";

    // Print the elements
    for (auto n : iota(0u, table.nrows())) {
        for (auto item : table.row(n))
            output << setw(table.get_width()) << item;
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
    std::cout << "Writing contents onto file '" << filename << "'\n";

    std::fstream file{string(filename), file.out};

    // Set precision
    file << std::setprecision(precision);

    auto inner_ids = column_ids | views::take(ncols()-1);
    auto last_id = column_ids.back();

    // Print headers
    for (auto id : inner_ids)
        file << id << ",";
    file << last_id << "\n";

    // Print columns
    for (auto n : iota(0u, this->nrows())) {
        auto row = this->row(n);
        auto inner_items = row | views::take(this->ncols() - 1);
        auto last_item = row.back();
        for (auto item : inner_items)
            file << item << ",";
        file << last_item << "\n";
    }
}

}

#endif
