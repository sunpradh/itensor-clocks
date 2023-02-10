#ifndef __UTILS_DATAFRAME_H
#define __UTILS_DATAFRAME_H

#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <concepts>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <functional>

using std::string;
using std::string_view;
using std::cout;
using std::setw;

template<typename T>
concept Container = requires(T a) {
    { std::begin(a) };
    { std::end(a) };
    { std::size(a) };
};

namespace utils {


template<Container T>
class Table {
    std::vector<string> column_names{};
    std::unordered_map<string, T> columns{};

    uint precision = 8;
    uint col_width = 12;
    uint n_rows = 0;

public:
    // Constructors
    Table(string_view col_name, const T & column);

    template<typename ... Args>
    Table(string_view first_col_name, const T & first_column, const Args & ... other_cols);

    // Add columns
    Table add_columns(string_view col_name, const T & column);

    template<typename ... Args>
    Table add_columns(string_view first_col_name, const T & first_column, const Args & ... other_cols);

    Table set_width(uint n) {
        col_width = n;
        return *this;
    }

    Table set_precision(uint n) {
        precision = n;
        return *this;
    }

    void print() ;
    void to_csv(string_view filename) ;

    T operator[](string_view col_name) {
        return columns[col_name];
    };
};


template<Container T>
Table<T>
Table<T>::add_columns(string_view col_name, const T & column){
    if (column.size() != n_rows)
        throw std::invalid_argument("Wrong size container");
    column_names.emplace_back(string(col_name));
    columns.insert(std::make_pair(col_name, column));
    return *this;
};

template<Container T>
template<typename ... Args>
Table<T>
Table<T>::add_columns(string_view first_col_name, const T & first_column, const Args & ... other_cols)
{
    add_columns(first_col_name, first_column);
    if (sizeof...(other_cols) < 2 || sizeof...(other_cols) % 2 != 0)
        throw std::invalid_argument("Invalid number of arguments");
    add_columns(other_cols...);
    return *this;
}


template<Container T>
Table<T>::Table(string_view col_name, const T & column) {
    n_rows = column.size();
    add_columns(col_name, column);
}

template<Container T>
template<typename ... Args>
Table<T>::Table(string_view first_col_name, const T & first_column, const Args & ... other_cols) {
    n_rows = first_column.size();
    add_columns(first_col_name, first_column, other_cols...);
}


template<Container T>
void
Table<T>::print() {
    // Set precision
    cout << std::setprecision(precision);

    // Print the headers
    for (auto col_name : column_names)
    cout << setw(col_width) << col_name;
    cout << "\n";

    // Print hrule
    for (uint i=0; i < col_width * column_names.size(); i++)
        cout << "-";
    cout << "\n";

    // Print the elements
    for (uint i=0; i < n_rows; i++) {
        for (auto col_name : column_names)
        cout << setw(col_width) << columns[col_name][i];
        cout << "\n";
    }
    cout << std::endl;
}


template<Container T>
void
Table<T>::to_csv(string_view filename) {
    std::fstream file{string(filename), file.out};

    // Set precision
    file << std::setprecision(precision);

    // Print headers
    for (auto name = column_names.begin(); name != column_names.end()-1; name++)
        file << *name << ",";
    file << column_names.back() << "\n";

    // Print columns
    for (uint i=0; i < n_rows; i++) {
        for (auto name = column_names.begin(); name != column_names.end()-1; name++)
            file << columns[*name][i] << ",";
        file << columns[column_names.back()][i] << "\n";
    }
}

}

#endif
