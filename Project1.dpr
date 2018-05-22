program Project1;

uses
  Vcl.Forms,
  ean_1 in 'ean_1.pas' {Form1},
  Bairstow in 'Bairstow.pas';

{$R *.res}

begin
  Application.Initialize;
  Application.MainFormOnTaskbar := True;
  Application.CreateForm(TForm1, Form1);
  Application.Run;
end.
