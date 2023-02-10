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
    // Constructors
    Table(string_view id, const T & column);

    template<typename ... Args>
    Table(string_view first_id, const T & first_column, const Args & ... other_cols);

    // Add columns
    Table add_columns(string_view id, const T & column);

    template<typename ... Args>
    Table add_columns(string_view first_id, const T & first_column, const Args & ... other_cols);

    Table set_width(unsigned int n) {
        col_width = n;
        return *this;
    }

    Table set_precision(unsigned int n) {
        precision = n;
        return *this;
    }

    auto ncols() { return columns.size(); }
    auto nrows() { return columns.begin()->second.size(); }
    auto size()  { return ncols() * nrows(); }

    void print() ;
    void to_csv(string_view filename) ;

    T operator[](string_view id) {
        return columns[id];
    };
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
Table<T>::print() {
    // Set precision
    cout << std::setprecision(precision);

    // Print the headers
    for (auto id : column_ids)
    cout << setw(col_width) << id;
    cout << "\n";

    // Print hrule
    for (auto i : iota(0u, col_width * ncols()))
        cout << "-";
    cout << "\n";

    // Print the elements
    for (auto i : iota(0u, nrows())) {
        for (auto id : column_ids)
        cout << setw(col_width) << columns[id][i];
        cout << "\n";
    }
    cout << std::endl;
}


template<ColumnType T>
void
Table<T>::to_csv(string_view filename) {
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
    for (auto i : iota(0u, nrows())) {
        for (auto id : inner_ids)
            file << columns[id][i] << ",";
        file << columns[last_id][i] << "\n";
    }
}

}

#endif
