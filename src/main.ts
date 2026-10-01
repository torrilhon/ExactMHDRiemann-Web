import "./ui/style.css";
import { downloadCsv, startApp } from "./ui/app.ts";

startApp();
document.getElementById("dl-csv")?.addEventListener("click", downloadCsv);
