"""`vhirprimerlab web`: serve the web app on this computer."""

from __future__ import annotations

import argparse
import threading
import webbrowser


def main(argv: list[str]) -> int:
    p = argparse.ArgumentParser(prog="vhirprimerlab web",
                                description="Run the vhirprimerlab web app on this computer.")
    p.add_argument("--host", default="127.0.0.1", help="address to listen on")
    p.add_argument("--port", type=int, default=8000, help="port to listen on")
    p.add_argument("--no-browser", action="store_true", help="do not open a browser window")
    a = p.parse_args(argv)

    import uvicorn

    url = f"http://{a.host}:{a.port}/"
    if not a.no_browser:
        threading.Timer(1.5, webbrowser.open, args=(url,)).start()
    print(f"vhirprimerlab web app on {url} (Ctrl+C to stop)")
    uvicorn.run("vhirprimerlab.web.app:app", host=a.host, port=a.port)
    return 0
