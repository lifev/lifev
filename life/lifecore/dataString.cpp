#include "dataString.hpp"

DataString::DataString(string str,int val,string help):
  _str(str),_val(val),_help(help)
{}

DataStringList::DataStringList(string title):
  _title(title)
{}

void DataStringList::add(string str,int val,string help)
{
  _list.push_back(DataString(str,val,help));
}

void DataStringList::showMe(ostream& c,bool val) const
{
  c << _title << " : " << endl;
  for(vector<DataString>::const_iterator ds = _list.begin();ds != _list.end();ds++){
    c << "   "  << ds->str() <<  " : " << ds->help();
    if(val) c << " (" << ds->val() << ")";
    c << endl;
  }  
}

int DataStringList::value(const string& str) const
{
  vector<DataString>::const_iterator ds = _list.begin();
  while(ds != _list.end()){
    if(ds->str() == str) return ds->val();
    ds++;
  };
  cout << "Error : " << str <<
    " is not in the list of possible choices for '" << _title << "'\n";
  showMe();
  exit(1);
  return 0;
}

