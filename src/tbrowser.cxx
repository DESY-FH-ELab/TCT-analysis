#include <TApplication.h>
#include <TBrowser.h>

int main(int argc, char ** argv)
{
   TApplication app("test_app", &argc, argv);
   TBrowser *b = new TBrowser("Results Browser");

   app.Run();
   return 0;
}

