/** these two functions (c-styled) are necessary for the 'sort' c/c++ library use */

/****************************************
   order the array by largest value 1st
   [first array element]
   proceeding down to smallest value
   [last array element]
****************************************/
bool orderByLargeToSmall(const DI &A, const DI& B)
{
   if (A.d < B.d)
      return (false);

   if (A.d > B.d)
      return (true);

   /** if A.d == B.d then ... */
   if (A.i < B.i)
      return (true);

   return (false);
}
/****************************************
   order the array by smallest value 1st
   [first array element]
   proceeding down to largest value
   [last array element]
****************************************/
bool orderBySmallToLarge(const DI &A, const DI& B)
{
   if (A.d < B.d)
      return (true);

   if (A.d > B.d)
      return (false);

   if (A.i < B.i)
      return (false);

   return (true);
}
