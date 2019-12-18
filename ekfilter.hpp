#define DEFINE_has_member(member_name)                                         \
    template <typename T>                                                      \
    class temp_has_member_##member_name                                        \
    {                                                                          \
        typedef char yes_type;                                                 \
        typedef long no_type;                                                  \
        template <typename U> static yes_type test(decltype(&U::member_name)); \
        template <typename U> static no_type  test(...);                       \
    public:                                                                    \
        static constexpr bool value = sizeof(test<T>(0)) == sizeof(yes_type);  \
    }

/// Shorthand for testing if "class_" has a member called "member_name"
///
/// @note "DEFINE_has_member(member_name)" must be used
///       before calling "has_member(class_, member_name)"
#define temp_has_member(class_, member_name)  temp_has_member_##member_name<class_>::value






#include<vector>
#include <iostream> // cout, endl
// #include <iomanip>  // std::boolalpha

typedef struct CubeSphereObject
{
    double x;
    double y;
    double z;
    double width;
    double length;
    double height;
    double heading;
} CubeSphereObject;

struct B
{
    bool heading;
};
struct C
{
    bool headings;
};
using std::cout;
using std::endl;
DEFINE_has_member(heading);

template<class T>
class test
{
public:
    test(){
        // check the existence of "sayHi"
        cout << "has_member(T, heading) " << temp_has_member(typename T:: value_type, heading) << endl;
        cout << endl;
    }
    ~test(){}
    
};

int main()
{
    // cout << std::boolalpha;  // display "true" or "false" for booleans
    test<std::vector<CubeSphereObject> > a;
    test<std::vector<B> > b;
    test<std::vector<C> > c;
    return 0;
}
